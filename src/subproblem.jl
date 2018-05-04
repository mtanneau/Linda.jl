
"""
    SubProblem is any optimization subproblem taking dual information π and σ
    and building new columns as a matrix and . 
"""
abstract type AbstractSubProblem end

"""
    solve for a subproblem, with both duals, should return:
    * costs returns a vector of N costs, one for each new column created
    * columns returns a MxN matrix (if the master problem has M constraints)
    * status indicates the status of the subproblem, must be an AbstractStatus
"""
function solve(::AbstractSubProblem, π::V1,σ::V2) where {V1<:AbstractVector{N1}, V2<:AbstractVector{N2}} where {N1<:Real, N2<:Real}
    costs = [0.0]
    columns = zeros(2,1)
    warn("Implement solve for the SubProblem")
    status = StatusError()
    return (status, costs, columns)
end
