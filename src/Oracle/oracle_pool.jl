"""
    LindaOraclePool

Pool of multiple sub-problems.


"""
mutable struct LindaOraclePool <: AbstractLindaOracle
    "Number of oracles in the pool."
    n::Int  # Number of oracles (i.e. sub-problems)

    new_columns::Set{Column}
    sp_dual_bound::Float64
    oracles::Vector{AbstractLindaOracle}

    LindaOraclePool(oracles) = new(size(oracles,1), Set{Column}(), -Inf, oracles)
end

"""
    Compute a set of columns with negative reduced cost.

# Arguments
- `pool`: The pool of oracles to be called
- `π`: Shadow prices (dual variables)
- `σ`: Shadow marginal cost (dual variables)
- `farkas`: Whether to perform Farkas pricing
- `tol_reduced_cost`: Numerical tolerance for reduced costs

# Returns
- `cols::Vector{Column}`
"""
function call_oracle!(
    pool::LindaOraclePool,
    π::AbstractVector{T1},
    σ::AbstractVector{T2};
    farkas::Bool=false,
    tol_reduced_cost::Float64=1.0e-6
) where{T1<:Real, T2<:Real}

    pool.new_columns = Set{Column}()

    # Check dimensions
    pool.n == size(σ, 1) || throw(DimensionMismatch(
        "Pool has $(pool.n) sub-problems but σ has size $(size(σ, 1))"
    ))

    pool.sp_dual_bound = (farkas ? (-Inf) : 0.0)
    best_red_cost = 0.0

    # price each sub-problem
    # go through sub-problems in random order
    nsp_priced = 0
    perm = randperm(pool.n)
    for r in perm

        o = pool.oracles[r]
        # solve sub-problem
        call_oracle!(
            o, π, σ[r],
            farkas=farkas,
            tol_reduced_cost=tol_reduced_cost
        )
        nsp_priced += 1

        # Check for infeasible sub-problem
        s = get_oracle_status(o)
        if s == PrimalInfeasible
            pool.status = PrimalInfeasible
            warn("Infeasible sub-problem")
            return nothing
        end

        # Update Lagrange bound
        pool.sp_dual_bound += get_sp_dual_bound(o)

        # Get new columns
        # TODO: parametrize
        cols = get_new_columns(o)
        for col in cols
            col.idx_subproblem = r
            rc = get_reduced_cost(col, π, σ[r], farkas)
            if rc < -tol_reduced_cost
                push!(pool.new_columns, col)
            end
            best_red_cost = min(rc, best_red_cost)

            if length(pool.new_columns) >= (0.1*pool.n) && (nsp_priced < 0.50*pool.n)
                # partial pricing
                # println("\tEarly pricing return")
                pool.sp_dual_bound = -Inf
                return nothing
            end
        end
        
    end

    # println("\tNCols: ", length(pool.new_columns))
    # println("\trc=: ", best_red_cost)
    # println("\tPriced:", nsp_priced)

    return nothing
     
end

get_new_columns(pool::LindaOraclePool) = pool.new_columns

get_sp_dual_bound(pool::LindaOraclePool) = pool.sp_dual_bound