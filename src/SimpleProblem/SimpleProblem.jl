"""
    SimpleProblem contains a full battery-included implementation of Linda Master and SubProblems
"""
module SimpleProblem

import MathProgBase

import Linda:
    AbstractSubProblem, AbstractMasterProblem, Column, PricingResult,
    solve_pricing,
    find_status, StatusError, StatusOptimal, StatusUnbounded, StatusInfeasible, ok, isinfeasible

export SimpleSubProblem, SimpleMasterProblem

#=======================
#
#   SimpleSubProblem
#
========================#

"""
    SimpleSubProblem
    Concrete implementation of AbstractSubProblem
"""
mutable struct SimpleSubProblem{N1<:Real,N2<:Real,N3<:Real,N4<:Real,N5<:Real} <: AbstractSubProblem
    
    # data for sub-problem
    costs::AbstractVector{N1}
    A::AbstractMatrix{N2}
    sense::AbstractVector{Char}
    b::AbstractVector{N3}
    vartypes::AbstractVector{Symbol}
    lb::AbstractVector{N4}
    ub::AbstractVector{N5}
    solver::MathProgBase.AbstractMathProgSolver
    columns::AbstractVector{Column}  # pool of columns
end

SimpleSubProblem(costs, A, sense, b, vartypes, lb, ub, solver) = SimpleSubProblem(costs, A, sense, b, vartypes, lb, ub, solver, Column[])

# TODO build convenient constructor functions for SimpleSubProblem

"""
    solve_pricing implementation for SimpleSubProblem
"""
function solve_pricing!(sp::SimpleSubProblem, π::V1, σ::V2, farkas_pricing=false) where {V1<:AbstractVector{N1}, V2<:AbstractVector{N2}} where {N1<:Real, N2<:Real}

    # dimension check
    size(π, 1) == size(sp.costs, 1) || DimensionMismatch("π's dimension is $(size(π, 1)) but should be $(size(sp.costs, 1))")

    # form perturbed objective for the sub-problem
    if farkas_pricing
        obj_ = - π
    else
        obj_ = sp.costs - π
    end

    # solve sub-problem
    result = MathProgBase.mixintprog(
        obj_,
        sp.A, sp.sense, sp.b, sp.vartypes, sp.lb, sp.ub, sp.solver
    )
    sp_status = find_status(result.status)

    # get solution
    if !ok(sp_status)
        # Error when solving sub-problem: return no column
        return PricinStatus(sp_status, Column[])
        # TODO: handle unbounded sub-problem
        #   In this case, a (primal) extreme ray is added to the master problem
    end
    
    # get column
    col = result.sol  # column
    cost = dot(sp.costs, col)  # native cost
    rc = result.objval - σ[1]  # reduced cost

    # check that reduced cost is indeed negative
    if rc < - 10.0^-6  # default solver tolerance for reduced costs
        columns = [Column(cost, col, true, false)]
    else
        columns = Column[]  # reduced cost is 0: return no column
    end

    return PricinStatus(sp_status, columns)
end

"""
    get_new_columns
    Return the set of columns that were generated during the last pricing.
"""
function get_new_columns(sp::SimpleSubProblem)
    return sp.columns
end


#=========================
#
#   SimpleMasterProblem
#
==========================#

"""
    SimpleMasterProblem
    Concrete implementation of MasterProblem.
"""
mutable struct SimpleMasterProblem{ST<:AbstractSubProblem} where {N1<:Number,N2<:Number}
    A::AbstractMatrix{N1}  # constraint matrix in original formulation
    b::AbstractVector{N2}  # right-hand side of linkin constraints

    rmp::MathProgBase.AbstractLinearQuadraticModel  # restricted Master Problem
    sp::ST  # sub-problem
end

# TODO build convenient constructor functions for SimpleMasterProblem
function SimpleMasterProblem(
    A::AbstractMatrix{N1},
    senses::AbstractVector{Char},
    b::AbstractVector{N2},
    solver::MathProgBase.AbstractMathProgSolver,
    sp::ST
) where {N1<:Number, N2<:Number, ST<:AbstractSubProblem}

    # dimension check
    nlinkingconstrs = size(A, 1)  # number of linking constraints
    nlinkingconstrs == size(senses, 1) || DimensionMismatch("")
    nlinkingconstrs == size(b, 1) || DimensionMismatch("")
    # instanciate rmp
    rmp = MathProgBase.LinearQuadraticModel(solver)

    # create empty constraints with right-hand sides
    for j=1:nlinkingconstrs
        if senses[j] == '='
            MathProgBase.addconstr!(rmp, [], [], b[j], b[j])
        elseif senses[j] == '<'
            MathProgBase.addconstr!(rmp, [], [], -Inf, b[j])
        elseif senses[j] == '>'
            MathProgBase.addconstr!(rmp, [], [], b[j], Inf)
        else
            error("Invalid input: senses[$(j)]=$(senses[j])")
        end
    end
    MathProgBase.addconstr!(rmp, [], [], 1.0, 1.0)  # convexity constraint

    return SimpleMasterProblem(A, b, rmp, sp)
end

"""
    subproblem
    returns the SubProblem object
"""
function subproblem(mp::SimpleMasterProblem{ST}) where {ST<:AbstractSubProblem}
    return mp.sp
end

"""
    compute_dual_variables!
    Compute the next dual iterate, by solving the RMP
"""
function compute_dual_variables!(mp::SimpleMasterProblem{ST}) where {ST<:AbstractSubProblem}

    # solve RMP
    MathProgBase.optimize!(mp.rmp)
    rmp_status = find_status(MathProgBase.status(mp.rmp))

    # check RMP status
    nlinkingconstrs = size(mp.A, 1)
    if !ok(rmp_status)
        # unexpected status when solving RMP
        # TODO: handle dual unboundedness in RMP
        π = zeros(0,)
        σ = zeros(1,)
    else
        # RMP solved to optimality
        dualsolution = MathProgBase.getconstrdual(mp.rmp)
        π = (mp.A)' * dualsolution[1:nlinkingconstrs]  # project dual vector to original space
        σ = dualsolution[nlinkingconstrs+1:end]
    end

    return (rmp_status, π, σ)
end

function add_columns!(mp::SimpleMasterProblem{ST}, columns::Vector{Column}) where ST<:AbstractSubProblem
    ncols = size(columns, 1)
    ncolsadded = 0

    constridx = collect(1:(1+size(mp.A, 1))  # assume only one sub-problem
    for column in columns
        # add column (assumed to have negative reduced cost)
        if column.isactive
            continue
        end
        ncolsadded += 1
        if column.isvertex
            # extreme vertex
            constrcoeff = vcat(mp.A * column.col, [1.0])
        else
            # extreme ray
            constrcoeff = vcat(mp.A * column.col, [0.0])
        end
        MathProgBase.addvar!(mp.rmp, constridx, constrcoeff, 0.0, Inf, col.cost)
        column.isactive = true
    end
    return ncolsadded
end

end  # module
