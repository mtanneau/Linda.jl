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
    MasterSolution
    Current primal-dual iterate in the master
"""
mutable struct MasterSolution

    # status of the Restricted Master Problem
    status::AbstractStatus

    # current dual iterate
    π::AbstractVector{Real}
    σ::AbstractVector{Real}

    # TODO: store primal iterate
    #   best is probably to store current basis

    # TODO: dual (Lagrange) bound, best upper bound (integer), best solution...

end

"""
    compute_dual_variables!
    
Computes the next dual iterate,
where `pi` is the dual variable associated to linking constraints,
and `sigma` is associated to the convexity constraint.
"""
function compute_dual_variables!(mp::AbstractMasterProblem)
    warn("Implement compute_dual_variables! for concrete MasterProblem types")
    rmp_status = StatusError()
    π = zeros(0,)
    σ = zeros(1,)
    return MasterSolution(rmp_status, π, σ)
end

"""
    subproblem 

Return the subproblem attached to a master.
To be defined in the implementation
"""
function subproblem(::AbstractMasterProblem) end

"""
    add_columns! is used to (validate and) add columns and
    
corresponding costs to the restricted master problem
"""
function add_columns!(::AbstractMasterProblem, columns::Vector{Column})
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

    # Display log
    println("  Iter     ncols")

    while newcols < maxcols

        ncgiter += 1
        # TODO: display relevant info, e.g.:
        # iter, number of columns, primal/dual bounds, etc...
        print(@sprintf("%6d", ncgiter))
        print(@sprintf("%10d", newcols))
        println()

        # I. Dual update
        rmp_sol = compute_dual_variables!(mp)
        if !ok(rmp_sol.status)
            # exit if problem encountered during dual update
            warn("RMP status $(rmp_sol.status) currently not handled, terminate")
            return rmp_sol.status
        end

        # II. Pricing step
        pricingresult = solve_pricing(sp, rmp_sol.π, rmp_sol.σ)
        # check pricing status
        if isinfeasible(pricingresult.status)
            # sub-problem is infeasible: problem is infeasible
            warn("Infeasible sub-problem: problem is infeasible")
            return StatusInfeasible()

        elseif !ok(pricingresult.status)
            # Early return caused by error when solving sub-problem
            # TODO: expand handling of return status
            warn("Pricing status $(res_pricing.sp_status) currently not handled, terminate")
            return pricingresult.status
        end

        # III. Update Master formulation
        if length(pricingresult.columns) == 0
            # no columns added: current solution is optimal
            # TODO: handle case where no columns are generated after heuristic solve
            return StatusOptimal()
        end

        ncolsadded = add_columns!(mp, pricingresult.columns)
        newcols += ncolsadded
    end
    return StatusTimeout()
end
