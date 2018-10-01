"""
    colgen()

Column-Generation algorithm.
"""
function solve_colgen!(
    env::LindaEnv,
    mp::LindaMaster{RMP},
    oracle::Oracle.AbstractLindaOracle
) where{RMP<:MPB.AbstractMathProgModel}

    # Pre-optimization stuff
    n_cg_iter::Int = 0
    num_bar_iter::Int = 0
    num_splx_iter::Int = 0
    time_start = time()
    time_mp_total::Float64 = 0.0
    time_sp_total::Float64 = 0.0
    time_cg_total::Float64 = 0.0

    if env[Val{:verbose}] == 1
        println(" Itn    Primal Obj      Dual Obj        NCols        (MP)     (SP)     (CG)")
    end

    # Main CG loop
    while n_cg_iter < env[Val{:num_cgiter_max}]
        # Solve RMP, update dual variables
        t0 = time()
        solve_rmp!(mp)
        time_mp_total += time() - t0
        num_bar_iter += MPB.getbarrieriter(mp.rmp)
        num_splx_iter += MPB.getsimplexiter(mp.rmp)

        if mp.rmp_status == Optimal
            farkas=false
        elseif mp.rmp_status == PrimalInfeasible
            farkas=true
        elseif mp.rmp_status == PrimalUnbounded
            println("Master Problem is unbounded.")
            return mp.mp_status
        else
            error("RMP status $(mp.rmp_status) not handled. Exiting.")
            return mp.mp_status
        end

        # Price
        t0 = time()
        Oracle.query!(
            oracle, mp.π, mp.σ,
            farkas=farkas,
            tol_reduced_cost=env.tol_reduced_cost.val,
            num_columns_max=env.num_columns_max.val
        )
        time_sp_total += time() - t0

        cols = Oracle.get_new_columns(oracle)
        lagrange_lb = (
            dot(mp.π, mp.rhs_constr_link)
            + Oracle.get_sp_dual_bound(oracle)
        )  # Compute Lagrange lower bound
        mp.dual_bound = lagrange_lb > mp.dual_bound ? lagrange_lb : mp.dual_bound

        # Log
        # Iteration count
        time_cg_total = time() - time_start
        if env[Val{:verbose}] == 1
            @printf("%4d", n_cg_iter)
            # Primal and Dual objectives
            @printf("%+18.7e", mp.primal_lp_bound)
            @printf("%+16.7e", mp.dual_bound)
            # RMP stats
            @printf("%10.0f", mp.num_columns_rmp)  # number of columns in RMP
            @printf("%9.2f", time_mp_total)
            @printf("%9.2f", time_sp_total)
            @printf("%9.2f", time_cg_total)
            @printf("%8d", MPB.getbarrieriter(mp.rmp))
            @printf("%6d", num_bar_iter)
            @printf("%10d", MPB.getsimplexiter(mp.rmp))
            @printf("%8d", num_splx_iter)
            print("\n")
        end

        # Check duality gap
        mp_gap = (
            abs(mp.primal_lp_bound - mp.dual_bound)
            / (1.0 + abs(mp.primal_lp_bound))
        )
        if mp_gap <= 10.0 ^-4
            mp.mp_status = Optimal
            
            time_cg_total = time() - time_start
            if env[Val{:verbose}] == 1
                println("Root relaxation solved.")
                println("Total time / MP: ", time_mp_total)
                println("Total time / SP: ", time_sp_total)
                println("Total time / CG: ", time_cg_total)
            end
            
            return mp.mp_status

        elseif farkas && length(cols) == 0
            
            mp.mp_status = PrimalInfeasible
            time_cg_total = time() - time_start
            if env[Val{:verbose}] == 1
                println("Master is infeasible.")
                println("Total time / MP: ", time_mp_total)
                println("Total time / SP: ", time_sp_total)
                println("Total time / CG: ", time_cg_total)
            end
            
            return mp.mp_status
        else
            # add columns
            add_columns!(mp, cols)
        end

        n_cg_iter += 1
    end

    time_cg_total = time() - time_start
    if env[Val{:verbose}] == 1
        println("Total time / MP: ", time_mp_total)
        println("Total time / SP: ", time_sp_total)
        println("Total time / CG: ", time_cg_total)
    end
    
    return mp.mp_status

end