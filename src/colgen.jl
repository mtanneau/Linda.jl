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
    if env[:verbose] == 1
        println(" Itn    Primal Obj      Dual Obj        NCols")
    end

    # Main CG loop
    while n_cg_iter < env[:num_cgiter_max]
        # Solve RMP, update dual variables
        solve_rmp!(mp)

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

        # Log
        # Iteration count
        if env[:verbose] == 1
            @printf("%4d", n_cg_iter)
            # Primal and Dual objectives
            @printf("%+18.7e", mp.primal_lp_bound)
            @printf("%+16.7e", mp.dual_bound)
            # RMP stats
            @printf("%10.0f", mp.num_columns_rmp)  # number of columns in RMP
            print("\n")
        end

        # Price
        Oracle.call_oracle!(env, oracle, mp.π, mp.σ, farkas=farkas)
        cols = Oracle.get_new_columns(oracle)
        lagrange_lb = (
            dot(mp.π, mp.rhs_constr_link)
            + Oracle.get_sp_dual_bound(oracle)
        )  # Compute Lagrange lower bound
        mp.dual_bound = lagrange_lb > mp.dual_bound ? lagrange_lb : mp.dual_bound

        # Check duality gap
        mp_gap = (
            abs(mp.primal_lp_bound - mp.dual_bound)
            / (1.0 + abs(mp.dual_bound))
        )
        if mp_gap <= 10.0 ^-4
            mp.mp_status = Optimal
            if env[:verbose] == 1
                println("Root relaxation solved.")
            end

            return mp.mp_status
        elseif farkas && length(cols) == 0
            if env[:verbose] == 1
                println("Master is infeasible.")
            end
            mp.mp_status = PrimalInfeasible
            return mp.mp_status
        else
            # add columns
            add_columns!(mp, cols)
        end

        n_cg_iter += 1
    end

    return mp.mp_status

end