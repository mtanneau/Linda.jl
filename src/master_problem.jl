"""
    MasterProblem is a main problem to solve using column generation
    Required methods to implement:  
    * `generate_column`
    * `subproblem()`
    * `add_column!`
    * `solve_restricted!`
    Optional methods (otherwise provided):
    * `solve!`
"""
abstract type AbstractMasterProblem{ST<:AbstractSubProblem} end

"""
    solve_restricted solves the MasterProblem with current columns
"""
function solve_restricted!(::AbstractMasterProblem)
    warn("Implement solve_restricted for concrete MasterProblem types")
    status = StatusError()
    π = [0.0]
    σ = [0.0]
    return (status, π, σ)
end


