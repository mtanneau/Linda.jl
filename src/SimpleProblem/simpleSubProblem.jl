"""
    SimpleSubProblem
A concrete implementation of AbstractSubProblem

# Attributes
-`costs::AbstractVector{Real}`: native objective in the sub-problem
-`A_link::AbstractMatrix{Real}`: matrix of the linking constraints,
    as given by the original formulation
-`A_sub::AbstractMatrix{Real}`: constraint matrix for the sub-problem
-`senses::Vector{Char}`: constraints' senses in the sub-problem
-`b::AbstractVector{Real}`: right-hand side of the sub-problem's constraints
-`vartypes::AbstractVector{Symbol}`: variables' types
-`lb::AbstractVector{Real}`: lower bounds on the original variables
-`ub::AbstractVector{Real}`: upper bounds on the original variables
-`solver::MathProgBase.AbstractMathProgSolver`: solver for the sub-problem
"""
mutable struct SimpleSubProblem{N1<:Real,N2<:Real,N3<:Real,N4<:Real,N5<:Real,N6<:Real} <: AbstractSubProblem
    
    # data for sub-problem
    costs::AbstractVector{N1}
    A_link::AbstractMatrix{N2}
    A_sub::AbstractMatrix{N3}
    senses::Vector{Char}
    b::AbstractVector{N4}
    vartypes::AbstractVector{Symbol}
    lb::AbstractVector{N5}
    ub::AbstractVector{N6}
    solver::MathProgBase.AbstractMathProgSolver

    problemidx::Integer # Sub-problem's number
end

num_subproblems(sp::SimpleSubProblem) = 1

# TODO build convenient constructor functions for SimpleSubProblem

"""
    solve_pricing implementation for SimpleSubProblem
"""
function solve_pricing(sp::SimpleSubProblem, π::AbstractVector{N1}, σ::AbstractVector{N2}, farkas_pricing=false) where {N1<:Real, N2<:Real}

    # dimension check
    size(π, 1) == size(sp.A_link, 1) || DimensionMismatch("π's dimension is $(size(π, 1)) but should be $(size(sp.A_link, 1))")

    # form perturbed objective for the sub-problem
    if farkas_pricing
        obj_ = - sp.A_link' * π
    else
        obj_ = sp.costs - sp.A_link' * π
    end

    # solve sub-problem
    result = MathProgBase.mixintprog(
        obj_,
        sp.A_sub, sp.senses, sp.b, sp.vartypes, sp.lb, sp.ub, sp.solver
    )
    sp_status = find_status(result.status)

    # get solution
    if !ok(sp_status)
        # Error when solving sub-problem: return no column
        return PricingResult(sp_status, Column[])
        # TODO: handle unbounded sub-problem
        #   In this case, a (primal) extreme ray is added to the master problem
    end
    
    # get column
    col = result.sol  # column
    cost = dot(sp.costs, col)  # native cost
    rc = result.objval - σ[1]  # reduced cost

    # check that reduced cost is indeed negative
    if rc < - 10.0^-6  # default solver tolerance for reduced costs
        columns = [Column(cost, sp.A_link * col, true, false, sp.problemidx)]
        # println("\tPricing: found column, rc = $(rc)")
        # println("\t\t", col)
    else
        columns = Column[]  # reduced cost is 0: return no column
    end

    return PricingResult(sp_status, columns)
end