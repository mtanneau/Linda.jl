
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
    * `costs` is a N-dimensional vector that contains the native costs of
        the generated columns.
    * `columns` is a MxN matrix, that contains the `N` generated columns.
    * `status` indicates the status of the subproblem, must be an AbstractStatus
"""
function solve_pricing(::AbstractSubProblem, π::V1,σ::V2, farkas_pricing = false) where {V1<:AbstractVector{N1}, V2<:AbstractVector{N2}} where {N1<:Real, N2<:Real}
    costs = [0.0]
    columns = zeros(2,1)
    warn("Implement solve_pricing for the SubProblem")
    status = StatusError()
    return (status, costs, columns)
end
