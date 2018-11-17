"""
    colgen()

Column-Generation algorithm.
"""
function solve_colgen!(
    env::LindaEnv,
    mp::LindaMaster{RMP},
    oracle::Oracle.AbstractLindaOracle;
    cg_log::Dict=Dict()
) where{RMP<:MPB.AbstractMathProgModel}

    # Pre-optimization stuff
    n_cg_iter::Int = 0
    time_start = time()
    cg_log[:time_mp_total] = 0.0
    cg_log[:time_sp_total] = 0.0
    cg_log[:time_cg_total] = 0.0

    cg_log[:num_iter_bar] = 0  # Total number of barrier inner iterations
    cg_log[:num_iter_spx] = 0  # Total number of simplex inner iterations

    cg_log[:nsp_priced] = 0  # Number of calls to sub-problems

    if env[Val{:verbose}] == 1
        println("Itn     Primal Obj        Dual Obj         NCols    MP(s)    SP(s)   Tot(s)  BarIter  SpxIter")
    end

    # Main CG loop
    while n_cg_iter < env[Val{:num_cgiter_max}]
        # Solve RMP, update dual variables
        t0 = time()
        solve_rmp!(mp)
        cg_log[:time_mp_total] += time() - t0

        if mp.rmp_status == Optimal
            farkas=false
        elseif mp.rmp_status == PrimalInfeasible
            farkas=true
        elseif mp.rmp_status == PrimalUnbounded
            env[Val{:verbose}] == 1 && println("Master Problem is unbounded.")
            break
        else
            @warn("RMP status $(mp.rmp_status) not handled.")
            break
        end

        # Price
        t0 = time()
        Oracle.query!(
            oracle, mp.π, mp.σ,
            farkas=farkas,
            tol_reduced_cost=env.tol_reduced_cost.val,
            num_columns_max=env.num_columns_max.val,
            log=cg_log
        )
        cg_log[:time_sp_total] += time() - t0

        cols = Oracle.get_new_columns(oracle)
        lagrange_lb = (
            dot(mp.π, mp.rhs_constr_link)
            + Oracle.get_sp_dual_bound(oracle)
        )  # Compute Lagrange lower bound
        mp.dual_bound = lagrange_lb > mp.dual_bound ? lagrange_lb : mp.dual_bound

        cg_log[:num_iter_bar] += MPB.getbarrieriter(mp.rmp)
        cg_log[:num_iter_spx] += MPB.getsimplexiter(mp.rmp)

        # Log
        # Iteration count
        if env[Val{:verbose}] == 1
            @printf("%4d", n_cg_iter)
            # Primal and Dual objectives
            @printf("%+18.7e", mp.primal_lp_bound)
            @printf("%+16.7e", mp.dual_bound)
            # RMP stats
            @printf("%10.0f", mp.num_columns_rmp)  # number of columns in RMP
            @printf("%9.2f", cg_log[:time_mp_total])
            @printf("%9.2f", cg_log[:time_sp_total])
            @printf("%9.2f", time() - time_start)
            @printf("%9d", cg_log[:num_iter_bar])
            @printf("%9d", cg_log[:num_iter_spx])
            print("\n")
            flush(Base.stdout)
        end

        # Check duality gap
        mp_gap = (
            abs(mp.primal_lp_bound - mp.dual_bound)
            / (1.0 + abs(mp.primal_lp_bound))
        )
        if mp_gap <= 10.0 ^-4
            mp.mp_status = Optimal
            
            # time_cg_total += time() - time_start
            if env[Val{:verbose}] == 1
                println("Root relaxation solved.")
            end
            
            break

        elseif farkas && length(cols) == 0
            
            mp.mp_status = PrimalInfeasible
            # time_cg_total += time() - time_start
            if env[Val{:verbose}] == 1
                println("Problem is infeasible.")
            end
            break
        else
            # add columns
            add_columns!(mp, cols)
        end

        n_cg_iter += 1
    end

    cg_log[:time_cg_total] = time() - time_start

    if env[Val{:verbose}] == 1
        @printf("Total time / MP: %.2fs\n", cg_log[:time_mp_total])
        @printf("Total time / SP: %.2fs\n", cg_log[:time_sp_total])
        @printf("Total time / CG: %.2fs\n", cg_log[:time_cg_total])
        @printf("Inner barrier iterations: %d\n", cg_log[:num_iter_bar])
        @printf("Inner simplex iterations: %d\n", cg_log[:num_iter_spx])
        @printf("Pricing calls: %d\n", cg_log[:nsp_priced])
    end
    
    return mp.mp_status

end