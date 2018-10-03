import Tulip
import Gurobi

function set_tolparam!(m::Tulip.TulipMathProgModel, ϵ)
    m.inner.env[:barrier_tol_conv] = ϵ
    return nothing
end
function set_tolparam!(m::Gurobi.GurobiMathProgModel, ϵ)
    Gurobi.setparam!(m.inner, "BarConvTol", ϵ)
    return nothing
end


"""
    pdcgm(env, mp, oracle)

Solve Master Problem using Primal-Dual Column Generation
"""
function pdcgm!(
    env::LindaEnv,
    mp::LindaMaster{RMP},
    oracle::Oracle.AbstractLindaOracle
) where{RMP<:MPB.AbstractMathProgModel}

    num_inneriter::Int = 0
    num_outeriter::Int = 0

    time_start::Float64 = time()
    time_cg = 0.0
    time_mp = 0.0
    time_sp = 0.0

    ϵ = 1e-2
    
    if env[Val{:verbose}] == 1
        println(" Itn    Primal Obj      Dual Obj        NCols        (MP)     (SP)     (CG)")
    end

    while num_outeriter < env[Val{:num_cgiter_max}]

        t0 = time()
        # Set RMP tolerances
        set_tolparam!(mp.rmp, ϵ)

        # Solve RMP
        solve_rmp!(mp)
        time_mp += time() - t0
        num_inneriter += MPB.getbarrieriter(mp.rmp)

        if (
            mp.rmp_status == Optimal
            || mp.rmp_status == DualFeasible
            || mp.rmp_status == PrimalDualFeasible
        )
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
        time_sp += time() - t0

        cols = Oracle.get_new_columns(oracle)
        lagrange_lb = (
            dot(mp.π, mp.rhs_constr_link)
            + Oracle.get_sp_dual_bound(oracle)
        )  # Compute Lagrange lower bound
        mp.dual_bound = lagrange_lb > mp.dual_bound ? lagrange_lb : mp.dual_bound

        # Log
        # Iteration count
        time_cg = time() - time_start
        if env[Val{:verbose}] == 1
            @printf("%4d", num_outeriter)
            # Primal and Dual objectives
            @printf("%+18.7e", mp.primal_lp_bound)
            @printf("%+16.7e", mp.dual_bound)
            # RMP stats
            @printf("%10.0f", mp.num_columns_rmp)  # number of columns in RMP
            @printf("%9.2f", time_mp)
            @printf("%9.2f", time_sp)
            @printf("%9.2f", time_cg)
            @printf("%8d", MPB.getbarrieriter(mp.rmp))
            @printf("%6d", num_inneriter)
            @printf("%9.2e", ϵ)
            print("\n")
        end

        # Check duality gap
        mp_gap = (
            abs(mp.primal_lp_bound - mp.dual_bound)
            / (1.0 + abs(mp.primal_lp_bound))
        )
        ϵ = min(1e-2, mp_gap / 10.0)

        if mp_gap <= 10.0 ^-4
            mp.mp_status = Optimal
            
            time_cg_total = time() - time_start
            if env[Val{:verbose}] == 1
                println("Root relaxation solved.")
                println("Total time / MP: ", time_mp)
                println("Total time / SP: ", time_sp)
                println("Total time / CG: ", time_cg)
            end
            
            return mp.mp_status

        elseif farkas && length(cols) == 0
            
            mp.mp_status = PrimalInfeasible
            time_cg = time() - time_start
            if env[Val{:verbose}] == 1
                println("Master is infeasible.")
                println("Total time / MP: ", time_mp)
                println("Total time / SP: ", time_sp)
                println("Total time / CG: ", time_cg)
            end
            
            return mp.mp_status
        else
            # add columns
            add_columns!(mp, cols)
        end

        num_outeriter += 1
    end

    return mp.mp_status
end