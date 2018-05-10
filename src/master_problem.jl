"""
    MasterProblem is a main problem to solve using column generation
    Required methods to implement:
    * `subproblem()`
    * `add_columns!`
    * `compute_dual_variables!``
    Optional methods (otherwise provided):
    * `solve!`
"""
abstract type AbstractMasterProblem{ST<:AbstractSubProblem} end

"""
    compute_dual_variables! computes the next dual iterate,
    where pi is the dual variable associated to linking constraints,
    and sigma is associated to the convexity constraint.
"""
function compute_dual_variables!(mp::AbstractMasterProblem)
    warn("Implement compute_dual_variables! for concrete MasterProblem types")
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
    add_columns! is used to (validate and) add columns and
    corresponding costs to the restricted master problem
"""
function add_columns!(::AbstractMasterProblem, costs::AbstractVector,columns::AbstractMatrix) end

"""
    solve! has a default version for any MasterProblem
    It adds column(s) to the restricted master while a new solution can be found
    in the subproblem. `maxcols` can be used to limit the number of
    new columns computed
"""
function solve!(mp::AbstractMasterProblem; maxcols::Integer = 5000)
    (status, π, σ) = compute_dual_variables!(mp)
    if !ok(status)
        return status
    end
    sp = subproblem(mp)
    newcols = 0
    while newcols < maxcols
        (status, costs, columns) = solve_pricing(sp, π, σ)
        if !ok(status)
            # not ok, no new negative cost column was found, early return 
            return status
        end
        add_columns!(mp, costs, columns)
        newcols += 1
        (status, π, σ) = compute_dual_variables!(mp)
        if !ok(status)
            return status
        end
    end
    return StatusTimeout()
end
