"""
    colgen()

Column-Generation algorithm.
"""
function solve_colgen!(
    mp::LindaMaster{RMP}
) where{RMP<:MPB.AbstractMathProgModel}

    # Pre-optimization stuff
    n_cg_iter = 0
    println(" Itn    Primal Obj      Dual Obj        NCols")
    # Main CG loop
    while n_cg_iter < 100
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
        @printf("%4d", n_cg_iter)
        # Primal and Dual objectives
        @printf("%+18.7e", mp.primal_lp_bound)
        @printf("%+16.7e", mp.dual_bound)
        # RMP stats
        @printf("%10.0f", mp.num_columns_rmp)  # number of columns in RMP
        print("\n")

        # Price
        Oracle.call_oracle!(mp.oracle, mp.π, mp.σ, farkas=farkas)
        cols = Oracle.get_new_columns(mp.oracle)
        lagrange_lb = (
            dot(mp.π, mp.rhs_constr_link)
            + Oracle.get_sp_dual_bound(mp.oracle)
        )  # Compute Lagrange lower bound
        mp.dual_bound = lagrange_lb > mp.dual_bound ? lagrange_lb : mp.dual_bound

        # Check duality gap
        mp_gap = (
            abs(mp.primal_lp_bound - mp.dual_bound)
            / (1.0 + abs(mp.dual_bound))
        )
        if mp_gap <= 10.0 ^-4
            mp.mp_status = Optimal
            println("Root relaxation solved.")
            return mp.mp_status
        else
            # add columns
            add_columns!(mp, cols)
        end

        n_cg_iter += 1
    end

    return mp.mp_status

end