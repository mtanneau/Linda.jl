"""
    MasterProblem is a main problem to solve using column generation
    Required methods to implement:  
    * `generate_column`
    * `subproblem()`
    * `add_column!`
    * `solve_restricted`
    Optional methods (otherwise provided):
    * `solve!`
"""
abstract type AbstractMasterProblem{ST<:AbstractSubProblem} end

"""
    solve_restricted solves the MasterProblem with current columns
"""
function solve_restricted(::AbstractMasterProblem)
    warn("Implement solve_restricted for concrete MasterProblem types")
    status = StatusError()
    π = [0.0]
    σ = [0.0]
    return (status, π, σ)
end

"""
    subproblem returns the subproblem attached to a master
    to define in the implementation
"""
function subproblem(::AbstractMasterProblem) end

"""
    solve! has a default version for any MasterProblem
    It adds column(s) to the master while a new solution can be found
    in the subproblem. `maxcols` can be used to limit the number of
    new columns computed
"""
function solve!(mp::AbstractMasterProblem; maxcols::Integer = 5000)
    (status, π, σ) = solve_restricted(mp)
    if !ok(status)
        return status
    end
    sp = subproblem(mp)
    newcols = 0
    while newcols < maxcols
        (status, costs, columns) = solve(sp, π, σ)
        if !ok(status)
            # not ok, no new negative cost column was found, early return 
            return status
        end
        add_column!(mp,costs,columns)
        newcols += 1
        (status, π, σ) = solve_restricted(mp)
        if !ok(status)
            return status
        end
    end
    return StatusTimeout()
end
