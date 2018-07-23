import MathProgBase
const MPB = MathProgBase

include("column.jl")

"""
    ProblemStatus

These status codes are based on those of MathOptInterface
- `PrimalFeasible`: A primal feasible solution has been found
- `DualFeasible` : A dual-feasible solution has been found
- `PrimalDualFeasible`: A primal- and a dual-feasible solutions have been found
- `Optimal`: A provably optimal solution has been found
- `PrimalInfeasible`: Primal is proved to be infeasible
- `PrimalUnbounded`: Problem is proved to be unbounded
- `Unknown`: No feasible solution nor proof of infeasibility yet
"""
@enum ProblemStatus Unknown PrimalFeasible DualFeasible PrimalDualFeasible Optimal PrimalInfeasible PrimalUnbounded

findStatus(::Symbol) = Unknown
findStatus(::Type{Val{:Infeasible}}) = PrimalInfeasible
findStatus(::Type{Val{:Unbounded}}) = PrimalUnbounded
findStatus(::Type{Val{:Optimal}}) = Optimal

"""
    LindaMaster

"""
mutable struct LindaMaster{RMP<:MPB.AbstractLinearQuadraticModel}

    #===========================================================================
        RMP data
    ===========================================================================#
    
    #=
        RMP Columns

        Variables 1 to `2*m`, in the RMP are artificial variables, with `m` the
            number of linking constraints.
        Columns are indexed from `2*m+1` and onwards.

        - `num_columns_rmp`: current number of columns in the RMP, excluding
            artificial variables and columns that have been generated, but are 
            not currently in the RMP.
        - `λ`: vector of primal variables, excluding artificial variables.
            It is updated each time the RMP is modified.
        
    =#
    num_columns_rmp::Int  # Number of columns currently in RMP
    λ::Vector{Float64}    # Vector of primal variables (columns coefficients)

    #=
        RMP Constraints

        The RMP has two types of constraints:
        - Convexity (1:`num_constr_cvxty`), one per sub-problem
        - Linking (`num_constr_cvxty+1`:`num_constr_cvxty+num_constr_link`)

        Dual variables:
        - `π`: vector of dual variables associated to linking constraints,
            a.k.a the vector of shadow prices.
            If the RMP is infeasible, π is a (dual) infeasibility ray.
        - `σ`: vector of dual variables associated to convexity constraints,
            a.k.a the vector of shadow marginal costs.
    =#
    num_constr_cvxty::Int   # Number of convexity constraints
    num_constr_link::Int  # Number of linking constraints
    rhs_constr_link::Vector{Float64}  # Right-hand side of linking constraints
    π::Vector{Float64}  # Shadow prices
    σ::Vector{Float64}  # Shadow marginal costs
    
    #=
        Restricted Master Problem (RMP)

        The RMP is solved at each iteration of the Column-Generation to 
            generate a new dual iterate.
    =#
    rmp::RMP  # Restricted Master Problem

    # Current status of the RMP
    # As long as the RMP is infeasbile, Farkas pricing will be used
    rmp_status::ProblemStatus  # Status of Restricted Master Problem


    #===========================================================================
        Column-Generation
    ===========================================================================#

    #=
        Master Problem Status

        The Master Problem is:
            - Feasible if a feasible solution has been found in the RMP.
            - Optimal if solved to optimality, i.e. there exists a primal
                feasible solution, and pricing fails to identify a column of
                negative reduced cost.
            - Infeasible if proven infeasible, i.e. RMP is infeasible and
                Farkas pricing fails to cut the infeasibility ray.
            - Unbounded if the RMP is unbounded
    =#
    mp_status::ProblemStatus  # Status of MasterProblem
    
    primal_lp_bound::Float64  # Primal bound for the linear (DW) relaxation
    primal_ip_bound::Float64  # Primal bound for the integer problem
    dual_bound::Float64   # Lagrange dual bound
    dual_bound_estimate::Float64  # Estimate of the Lagrange dual bound

    column_pool::Vector{Column}  # Pool of all the generated columns

    # Create an empty Master
    function LindaMaster(
        rmp::RMP,
        num_constr_cvxty::Int,
        num_constr_link::Int,
        rhs_constr_link::AbstractVector{Float64},
        column_pool::Vector{Column}
    ) where RMP<:MPB.AbstractLinearQuadraticModel

        # Dimension check
        n = MPB.numvar(rmp)
        n == (2 * num_constr_link) || throw(ErrorException(
            "RMP has $(n) variables instead of $(2*num_constr_link)"
        ))
        m = MPB.numconstr(rmp)
        m == (num_constr_cvxty + num_constr_link) || throw(
            ErrorException(
                "RMP has $(m) constraints instead of $(num_constr_cvxty + num_constr_link)"
            )
        )
        num_constr_link == size(rhs_constr_link, 1) || throw(DimensionMismatch(
            "RMP has $(num_constr_link) linking constraints but b has size $(size(rhs_constr_link))"
        ))

        # Instanciate model
        mp = new{RMP}()

        mp.num_columns_rmp = 0
        mp.λ = Vector{Float64}(0)

        mp.num_constr_cvxty = num_constr_cvxty
        mp.num_constr_link = num_constr_link
        mp.rhs_constr_link = rhs_constr_link
        mp.π = Vector{Float64}(num_constr_link)
        mp.σ = Vector{Float64}(num_constr_cvxty)

        mp.rmp = rmp
        mp.rmp_status = Unknown

        mp.mp_status = Unknown
        mp.primal_lp_bound = Inf
        mp.primal_ip_bound = Inf
        mp.dual_bound = -Inf
        mp.dual_bound_estimate = -Inf

        # Initial pool of columns. These are not added to the RMP yet.
        mp.column_pool = column_pool

        return mp
    end
end

function solve!(master::LindaMaster; verbose=1)

    # Run-time initialization


    # main CG loop
    while master.num_columns_rmp < 10^4

        # Solve Restricted Master Problem to update dual variables
        solve_rmp!(master)

        # Log
        if verbose == 1
            println()
        end

        # Look for early termination
        if master.mp_status == PrimalUnbounded
            println("Master Problem is unbounded.")
            return nothing
        end

        # Pricing step


        
    end



    return nothing
end

"""
    solve_rmp!(master)

Solve Restricted Master Problem, and update current dual iterate.
"""
function solve_rmp!(master::LindaMaster)

    MPB.optimize!(master.rmp)
    rmp_status = findStatus(Val{MPB.status(master.rmp)})
    master.rmp_status = rmp_status

    # Update dual iterate
    if rmp_status == Optimal || rmp_status == PrimalDualFeasible
        # update primal bound
        master.primal_lp_bound = MPB.getobjval(master.rmp)
        master.mp_status = PrimalFeasible

        # update current column coefficients
        x = MPB.getsolution(master.rmp)
        master.λ .= x[(2*master.num_constr_link+1):end]

        # update dual variables
        y = MPB.getconstrduals(master.rmp)
        master.σ .= y[1:master.num_constr_cvxty]
        master.π .= y[(master.num_constr_cvxty+1):(master.num_constr_cvxty+master.num_constr_link)]
    
    elseif rmp_status == :PrimalInfeasible
        # update primal bound
        master.primal_lp_bound = Inf

        # update dual variables
        y = MPB.getinfeasibilityray(master.rmp)
        master.σ .= 0.0
        master.π .= y[(master.num_constr_cvxty+1):(master.num_constr_cvxty+master.num_constr_link)]

    elseif rmp_status == :PrimalUnbounded
        master.primal_bound = -Inf
        master.dual_bound = -Inf

        # RMP is unbounded, thus Master is unbounded
        master.mp_status = PrimalUnbounded

    else
        # TODO: raise error
        throw(ErrorException("Unhandled status: $(rmp_status)"))
    end

    return nothing
end