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
    π = zeros(0,)
    σ = zeros(1,)
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
function add_columns!(::AbstractMasterProblem, columns::Array{Column, 1})
    warn("Implement add_columns! for concrete MasterProblem types")
    return 0
end

"""
    solve! has a default version for any MasterProblem
    It adds column(s) to the restricted master while a new solution can be found
    in the subproblem. `maxcols` can be used to limit the number of
    new columns computed
"""
function solve!(mp::AbstractMasterProblem; maxcols::Integer = 5000)

    # TODO: initialize restricted master (e.g. artificial variables)
    # TODO: heuristic hotstart

    sp = subproblem(mp)
    newcols = 0  # number of columns added to the master problem
    ncgiter = 0  # number of Column Generation iterations

    while newcols < maxcols

        ncgiter += 1
        # TODO: display relevant info, e.g.:
        # iter, number of columns, primal/dual bounds, etc...

        # I. Dual update
        (rmp_status, π, σ) = compute_dual_variables!(mp)
        if !ok(rmp_status)
            # exit if problem encountered during dual update
            warn("RMP status $(rmp_status) currently not handled, terminate")
            return rmp_status
        end

        # II. Pricing step
        sp_status = solve_pricing(sp, π, σ)
        # check pricing status
        if isinfeasible(sp_status)
            # sub-problem is infeasible: problem is infeasible
            warn("Infeasible sub-problem: problem is infeasible")
            return StatusInfeasible()

        elseif !ok(sp_status)
            # Early return caused by error when solving sub-problem
            # TODO: expand handling of return status
            warn("Pricing status $(sp_status) currently not handled, terminate")
            return sp_status
        end

        # III. Update Master formulation
        columns = get_new_columns(sp)  # get columns from sub-problem
        if length(columns) == 0
            # no columns added: current solution is optimal
            return StatusOptimal()
        end

        ncolsadded = add_columns!(mp, columns)
        newcols += ncolsadded
    end
    return StatusTimeout()
end
