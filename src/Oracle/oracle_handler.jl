TOL_REDCOST = 10.0^-6  # numerical tolerance for reduced costs

"""

"""
mutable struct LindaOracleHandler <: AbstractLindaOracle
    n::Int  # Number of oracles (i.e. sub-problems)

    new_columns::Set{Column}
    sp_dual_bound::Float64
    oracles::Vector{AbstractLindaOracle}

    LindaOracleHandler(oracles) = new(size(oracles,1), Set{Column}(), -Inf, oracles)
end

function call_oracle!(
    env::LindaEnv,
    handler::LindaOracleHandler,
    π::AbstractVector{T1},
    σ::AbstractVector{T2};
    farkas::Bool=false
) where{T1<:Real, T2<:Real}
    # println("\tPricing, Farkas=$farkas")
    handler.new_columns = Set{Column}()

    # Check dimensions
    handler.n == size(σ, 1) || throw(DimensionMismatch(
        "Handler has $(handler.n) sub-problems but σ has size $(size(σ, 1))"
    ))

    handler.sp_dual_bound = (farkas ? (-Inf) : 0.0)
    best_red_cost = 0.0

    # price each sub-problem
    for (r, o) in enumerate(handler.oracles)
        # solve sub-problem
        call_oracle!(env, o, π, σ[r])

        # Check for infeasible sub-problem
        s = get_oracle_status(o)
        if s == PrimalInfeasible
            handler.status = PrimalInfeasible
            return nothing
        end

        # Update Lagrange bound
        handler.sp_dual_bound += get_sp_dual_bound(o)

        # Get new columns
        cols = get_new_columns(o)
        for col in cols
            col.idx_subproblem = r
            rc = get_reduced_cost(col, π, σ[r], farkas)
            if rc < -env[:tol_reduced_cost]
                push!(handler.new_columns, col)
            end
            best_red_cost = min(rc, best_red_cost)
        end
        
    end

    # println("\tNCols: ", length(handler.new_columns))
    # println("\trc=: ", best_red_cost)

    return nothing
     
end

get_new_columns(handler::LindaOracleHandler) = handler.new_columns

get_sp_dual_bound(handler::LindaOracleHandler) = handler.sp_dual_bound