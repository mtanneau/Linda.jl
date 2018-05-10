
"""
    SubProblem is any optimization subproblem taking dual information π and σ
    and building new columns as a matrix and . 
"""
abstract type AbstractSubProblem end

"""
    solve_pricing performs the pricing of the current dual iterate, and returns a set of `N` columns.
    Note that `N` may be equal to zero.
    
    Arguments:
    * `pi` is the vector of dual variables associated to linking constraints in the master
    * `sigma` is the (Vector of) dual variable(s) associated to convexity constraint(s)
    * `farkas_pricing` indicates whether to perform Farkas pricing, or regular pricing
    
    Returns:
    * `costs` is N-dimensional vector that contains the cost if each generated column
    * `columns` is a MxN matrix, where `M` is the number of linking constraints in the master,
    * `status` indicates the status of the subproblem, must be an AbstractStatus
"""
function solve_pricing(::AbstractSubProblem, π::V1,σ::V2, farkas_pricing = false) where {V1<:AbstractVector{N1}, V2<:AbstractVector{N2}} where {N1<:Real, N2<:Real}
    costs = [0.0]
    columns = zeros(2,1)
    warn("Implement solve_pricing for the SubProblem")
    status = StatusError()
    return (status, costs, columns)
end
