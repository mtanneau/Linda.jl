"""
    LindaOraclePool

Pool of multiple sub-problems.

# Attributes
- `n`: Number of oracles (i.e., sub-problems) in the pool
- `new_columns`: Set of columns (with negative reduced cost) that were generated
    by the last query.
- `status`: Pricing status of the current query
- `sp_dual_bound`: The best known dual bound on the optimal value of the sub-problem
- `oracles`: The list of oracles in the pool
"""
mutable struct LindaOraclePool <: AbstractLindaOracle
    n::Int  # Number of oracles (i.e. sub-problems)

    new_columns::Set{Column}
    status::Status
    sp_dual_bound::Float64

    oracles::Vector{AbstractLindaOracle}

    LindaOraclePool(oracles) = new(size(oracles,1), Set{Column}(), Unknown, -Inf, oracles)
end

"""
    query!(pool, π, σ; kwargs...)    
Compute a set of columns with negative reduced cost.

# Arguments
- `pool`: The pool of oracles to be called
- `π`: Shadow prices (dual variables)
- `σ`: Shadow marginal cost (dual variables)
- `farkas`: Whether to perform Farkas pricing
- `tol_reduced_cost`: Numerical tolerance for reduced costs
- `num_columns_max`: Maximum number of columns to be generated before pricing stops.
- `prop_sp_priced_max`: Maximum proportion of sub-problems to be solved before pricing can stop.

# Returns
- `status`: Termination status
"""
function query!(
    pool::LindaOraclePool,
    π::AbstractVector{T1},
    σ::AbstractVector{T2};
    farkas::Bool=false,
    tol_reduced_cost::Float64=1.0e-6,
    num_columns_max::Int=typemax(Int64),
    prop_sp_priced_max::Float64=1.0
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
        query!(
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
            return pool.status
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

            if length(pool.new_columns) >= num_columns_max
                # early return
                pool.sp_dual_bound = -Inf
                return pool.status
            end
        end
        
    end

    # Pricing ended because all sub-problems were solved to optimality
    pool.status = Optimal

    # println("\tNCols: ", length(pool.new_columns))
    # println("\trc=: ", best_red_cost)
    # println("\tPriced:", nsp_priced)

    return pool.status
     
end

get_new_columns(pool::LindaOraclePool) = pool.new_columns

get_sp_dual_bound(pool::LindaOraclePool) = pool.sp_dual_bound