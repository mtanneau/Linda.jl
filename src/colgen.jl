using LinearAlgebra

"""
    colgen()

Column-Generation algorithm.
"""
function solve_colgen!(
    env::LindaEnv,
    mp::Master,
    oracles::Vector{To};
    cg_log::Dict=Dict()
) where{To<:Oracle}

    # Pre-optimization stuff
    n_cg_iter::Int = 0
    time_start = time()
    cg_log[:time_mp_total] = 0.0
    cg_log[:time_sp_total] = 0.0
    cg_log[:time_cg_total] = 0.0

    cg_log[:num_iter_bar] = 0  # Total number of barrier inner iterations
    cg_log[:num_iter_spx] = 0  # Total number of simplex inner iterations

    cg_log[:nsp_priced] = 0  # Number of calls to sub-problems

    n_slow_progress = 0  # Number of consecutive iterations with slow progress

    if env[Val{:verbose}] == 1
        @printf "%4s  %15s  %15s  %6s  %8s  %8s  %8s  %8s\n" "Iter" "Primal Obj." " Dual bound" "%Gap" "#Cols" "Time(s)" "BarIter"  "SpxIter"
    end

    # CG loop
    while true

        #
        # Solve RMP
        # 
        cg_log[:time_mp_total] += @elapsed MOI.optimize!(mp.rmp)
        st = MOI.get(mp.rmp, MOI.TerminationStatus())

        cg_log[:num_iter_bar] += MOI.get(mp.rmp, MOI.BarrierIterations())
        cg_log[:num_iter_spx] += MOI.get(mp.rmp, MOI.SimplexIterations())

        # Check RMP termination status
        if st == MOI.OPTIMAL
            # RMP is primal-feasible
            z_old = mp.primal_lp_bound
            z_new = MOI.get(mp.rmp, MOI.ObjectiveValue())
            mp.primal_lp_bound = z_new

            # Optimality gap
            g = (
                abs(z_new - mp.dual_bound)
                / (1e-8 + abs(z_new))
            )
            isfinite(g) || (g = Inf)
            g <= 1e-4 && (mp.mp_status = MOI.OPTIMAL)
            mp.mp_gap = g

            # Check solution improvement
            if isfinite(z_new) && isfinite(z_old)
                z = abs(z_new - z_old) / (1e-8 + abs(z_old))
                # Relative improvement should be at least 0.01%
                if z > 1e-4 
                    n_slow_progress = 0
                else
                    n_slow_progress += 1
                end
            end

            # Use regular pricing
            farkas = false

        elseif st == MOI.INFEASIBLE
            # RMP is infeasible
            mp.primal_lp_bound = Inf

            # Use Farkas pricing
            farkas = true

        elseif st == MOI.DUAL_INFEASIBLE
            # RMP is unbounded, thus MP is unbounded
            mp.primal_lp_bound = -Inf
            mp.mp_status = MOI.DUAL_INFEASIBLE

        else
            error("RMP solver exited with status $st")
        end

        #
        # Log
        # 
        if env[Val{:verbose}] == 1
            @printf "%4d" n_cg_iter
            # Primal and Dual objectives
            @printf("  %+15.8e", mp.primal_lp_bound)
            @printf("  %+15.8e", mp.dual_bound)
            @printf("  %6.2f", mp.mp_gap)
            # RMP stats
            @printf("  %8d", length(mp.columns))  # number of columns in RMP
            # @printf("%9.2f", cg_log[:time_mp_total])
            # @printf("%9.2f", cg_log[:time_sp_total])
            @printf("  %8.2f", time() - time_start)
            @printf("  %8d", cg_log[:num_iter_bar])
            @printf("  %8d", cg_log[:num_iter_spx])
            @printf("\n")
        end

        # 
        # Check stopping criterion
        # 
        if mp.mp_status == MOI.DUAL_INFEASIBLE || mp.mp_status == MOI.OPTIMAL
            # stop
            break
        elseif n_cg_iter > env[Val{:num_cgiter_max}]
            # Maximum number of CG iteration reached
            mp.mp_status = MOI.ITERATION_LIMIT
            break
        elseif time() - time_start > env[:time_limit]
            # Time limit reached
            mp.mp_status = MOI.TIME_LIMIT
            break
        elseif length(mp.columns) > env.num_columns_max.val
            mp.mp_status = MOI.OTHER_LIMIT
            break
        elseif n_slow_progress > 10
            mp.mp_status = MOI.SLOW_PROGRESS
            break
        end

        # 
        # Pricing step
        # 
        # Check that dual status is OK
        dst = MOI.get(mp.rmp, MOI.DualStatus())
        if !(dst == MOI.FEASIBLE_POINT || dst == MOI.INFEASIBILITY_CERTIFICATE)
            env[Val{:verbose}] == 1 && println("Dual status in RMP is $(dst)")
            mp.mp_status = MOI.OTHER_ERROR
            break
        end

        # Recover dual variables
        for (i, cidx) in enumerate(mp.con_link)
            mp.π[i] = MOI.get(mp.rmp, MOI.ConstraintDual(), cidx)
        end
        for (i, cidx) in enumerate(mp.con_cvx)
            mp.σ[i] = MOI.get(mp.rmp, MOI.ConstraintDual(), cidx)
        end
        
        # TODO: column pool

        # Generate new columns
        new_cols = Tuple{Column, Float64}[]
        z_lagrange = farkas ? -Inf : dot(mp.π, mp.rhs_link)
        for (i, o) in enumerate(oracles)
            # TODO: Check stopping criterion
            
            # Update 
            update!(o, farkas, mp.π, mp.σ[i])

            # Solve sub-problem
            optimize!(o)

            # Recover columns
            cols = get_columns(o)
            append!(new_cols, cols)

            # Update Lagrange bound
            z_lagrange += get_dual_bound(o)
        end

        # Lagrange bound update
        mp.dual_bound = max(mp.dual_bound, z_lagrange)

        # Add new columns
        ncols_added = 0
        for (col, rc) in new_cols
            if rc <= - env.tol_reduced_cost.val
                add_column!(mp, col)
                ncols_added += 1
            end
        end
        
        if ncols_added == 0
            if st == MOI.OPTIMAL
                # Master is optimal
                mp.mp_status = MOI.OPTIMAL
            elseif st == MOI.INFEASIBLE
                # Master is infeasible
                mp.mp_status = MOI.INFEASIBLE
            end
            break
        end

        # Bump iteration count
        n_cg_iter += 1

    end  # CG loop

    # Logs
    cg_log[:time_cg_total] = time() - time_start
    cg_log[:dual_bound] = mp.dual_bound
    cg_log[:primal_bound] = mp.primal_lp_bound
    cg_log[:n_cg_iter] = n_cg_iter
    cg_log[:status] = mp.mp_status

    if env[Val{:verbose}] == 1
        println("CG terminated with status $(mp.mp_status)")
        println()
        @printf("Total time / MP: %.2fs\n", cg_log[:time_mp_total])
        @printf("Total time / SP: %.2fs\n", cg_log[:time_sp_total])
        @printf("Total time / CG: %.2fs\n", cg_log[:time_cg_total])
        @printf("Inner barrier iterations: %d\n", cg_log[:num_iter_bar])
        @printf("Inner simplex iterations: %d\n", cg_log[:num_iter_spx])
        @printf("Pricing calls: %d\n", cg_log[:nsp_priced])
    end
    
    return mp.mp_status

end