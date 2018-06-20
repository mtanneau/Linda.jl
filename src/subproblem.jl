
"""
    AbstractSubProblem 
    
Any optimization subproblem taking dual information π and σ to generate new
columns.
"""
abstract type AbstractSubProblem end

"""
    getprobindex(::AbstractSubProblem)

Return index of sub-problem
"""
function getprobindex(::AbstractSubProblem)
    warn("Implement getprobindex for concrete implementations")
    return 1
end


"""
    PricingResult 

Holds relevant information following a pricing step, including the return status
of the sub-problem, and the columns that were generated (if any).
"""
struct PricingResult
    status::AbstractStatus   # return status of the sub-problem
    columns::Vector{Column}  # columns generated during pricing
end

"""
    solve_pricing 

Perform pricing for the current dual iterate, and return a set of `N` columns,
where `N` may be equal to zero.
    
# Arguments:
-`::AbstractSubProblem`: Sub-problem
-`pi::AbstractVector{Real}`: Vector of dual variables (shadow prices)
-`sigma::Real`: Dual variable (shadow marginal cost)
-`farkas_pricing::Bool`: Whether to perform Farkas or regular pricing

# Returns:
-`result::PricingSTatus`: Optimization status for the sub-problem, and set of
    newly generated columns (if any)
"""
function solve_pricing(
    ::AbstractSubProblem,
    π::V1,
    σ::Real;
    farkas_pricing=false
) where {V1<:AbstractVector{N1}} where {N1<:Real}
    
    warn("Implement solve_pricing for concrete SubProblem types")
    status = StatusError()
    columns = Column[]  # empty set of columns

    return PricingResult(status, columns)
end
