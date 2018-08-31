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
    n_cg_iter = 0
    time_mp_total = 0.0
    time_sp_total = 0.0
    time_cg_total = 0.0

    tic()  # start

    if env[:verbose] == 1
        println(" Itn    Primal Obj      Dual Obj        NCols    Time (s)")
    end

    # Main CG loop
    while n_cg_iter < env[:num_cgiter_max]
        # Solve RMP, update dual variables
        tic()
        solve_rmp!(mp)
        time_mp_total += toq()

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
        tic()
        Oracle.call_oracle!(
            oracle, mp.π, mp.σ,
            farkas=farkas,
            tol_reduced_cost=env[:tol_reduced_cost]
        )
        time_sp_total += toq()

        cols = Oracle.get_new_columns(oracle)
        lagrange_lb = (
            dot(mp.π, mp.rhs_constr_link)
            + Oracle.get_sp_dual_bound(oracle)
        )  # Compute Lagrange lower bound
        mp.dual_bound = lagrange_lb > mp.dual_bound ? lagrange_lb : mp.dual_bound

        # Log
        # Iteration count
        if env[:verbose] == 1
            @printf("%4d", n_cg_iter)
            # Primal and Dual objectives
            @printf("%+18.7e", mp.primal_lp_bound)
            @printf("%+16.7e", mp.dual_bound)
            # RMP stats
            @printf("%10.0f", mp.num_columns_rmp)  # number of columns in RMP
            @printf("%9.2f", time_mp_total + time_sp_total)
            print("\n")
        end

        # Check duality gap
        mp_gap = (
            abs(mp.primal_lp_bound - mp.dual_bound)
            / (1.0 + abs(mp.primal_lp_bound))
        )
        if mp_gap <= 10.0 ^-4
            mp.mp_status = Optimal
            
            time_cg_total += toq()
            if env[:verbose] == 1
                println("Root relaxation solved.")
                println("Total time / MP: ", time_mp_total)
                println("Total time / SP: ", time_sp_total)
                println("Total time / CG: ", time_cg_total)
            end
            
            return mp.mp_status

        elseif farkas && length(cols) == 0
            
            mp.mp_status = PrimalInfeasible
            time_cg_total += toq()
            if env[:verbose] == 1
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

    time_cg_total += toq()
    if env[:verbose] == 1
        println("Total time / MP: ", time_mp_total)
        println("Total time / SP: ", time_sp_total)
        println("Total time / CG: ", time_cg_total)
    end
    
    return mp.mp_status

end