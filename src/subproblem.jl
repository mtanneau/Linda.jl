
"""
    SubProblem is any optimization subproblem taking dual information π and σ
    and building new columns as a matrix and . 
"""
abstract type AbstractSubProblem end

"""
    price performs the pricing of the current dual iterate, and returns a set
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
    return status
end

"""
    get_new_columns
    Returns a (possibly empty) array of Column objects
"""
function get_new_columns(sp::AbstractSubProblem)

    warn("Implement get_new_columns for concrete SubProblem types")
    columns = Array{Column, 1}()  # empty set of columns
    return columns

end
