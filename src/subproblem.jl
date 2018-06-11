
"""
    SubProblem is any optimization subproblem taking dual information π and σ
    and building new columns as a matrix and . 
"""
abstract type AbstractSubProblem end

"""
    PricingResult stores relevant information following a pricing step,
    including the return status of the sub-problem, and the columns that were
    generated (if any).
"""
struct PricingResult
    status::AbstractStatus  # return status of the sub-problem
    columns::Vector{Column}  # columns generated during pricing
end

"""
    num_subproblems

Return number of (independant) sub-problems.
"""
function num_subproblems(::AbstractSubProblem)
    warn("Implement num_subproblems for concrete SubProblem types")
    return 0
end

"""
    solve_pricing performs the pricing of the current dual iterate, and returns a set
    of `N` columns, where `N` may be equal to zero.
    
    Arguments:
    * `pi` is the vector of dual variables associated to linking constraints
    * `sigma` is the (vector of) dual variable(s) associated to convexity constraint(s)
    * `farkas_pricing` indicates whether to perform Farkas or regular pricing
    
    Returns:
    * `status` indicates the status of the subproblem, must be an AbstractStatus
"""
function solve_pricing(::AbstractSubProblem, π::V1,σ::V2, farkas_pricing = false) where {V1<:AbstractVector{N1}, V2<:AbstractVector{N2}} where {N1<:Real, N2<:Real}
    
    warn("Implement solve_pricing for concrete SubProblem types")
    status = StatusError()
    columns = Column[]  # empty set of columns

    return PricingResult(status, columns)
end
