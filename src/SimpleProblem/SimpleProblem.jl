"""
    SimpleProblem contains a full battery-included implementation of Linda Master and SubProblems
"""
module SimpleProblem

import MathProgBase

import Linda:
    AbstractSubProblem, AbstractMasterProblem, solve_pricing,
    find_status, StatusError, StatusOptimal, StatusUnbounded, StatusInfeasible, ok

export SimpleSubProblem, SimpleMasterProblem

mutable struct SimpleSubProblem{N1<:Real,N2<:Real,N3<:Real,N4<:Real,N5<:Real} <: AbstractSubProblem
    costs::AbstractVector{N1}
    A::AbstractMatrix{N2}
    sense::AbstractVector{Char}
    b::AbstractVector{N3}
    vartypes::AbstractVector{Symbol}
    lb::AbstractVector{N4}
    ub::AbstractVector{N5}
    solver::MathProgBase.AbstractMathProgSolver
end

# TODO build convenient constructor functions for SimpleSubProblem

"""
    solve_pricing implementation for SimpleSubProblem
"""
function solve_pricing(sp::SimpleSubProblem, π::V1, σ::N2, farkas_pricing=false) where {V1<:AbstractVector{N1}} where {N1<:Real, N2<:Real}

    # form perturbed objective for the sub-problem
    obj = farkas_pricing ? (-π) : (sp.costs - π)

    # solve sub-problem
    result = MathProgBase.mixintprog(
        obj,
        sp.A, sp.sense, sp.b, sp.vartypes, sp.lb, sp.ub, sp.solver
    )
    final_status = find_status(result.status)

    # get solution
    if !ok(final_status)
        # Error when solving sub-problem: return no column
        return (final_status, zeros(0,), zeros(0,0))
        # TODO: handle unbounded sub-problem
        #   In this case, a (primal) extreme ray is added to the master problem
    end
    col = result.sol  # column
    cost = dot(sp.costs, col)  # native cost
    rc = result.objval - σ  # reduced cost

    # check that reduced cost is indeed negative
    if rc < -eps(Float32)
        return (final_status, [cost], hcat(col))
    else
        # reduced cost is 0: return no column
        return (final_status, zeros(0,), zeros(0,0))
    end
end

mutable struct SimpleMasterProblem{SimpleSubProblem,N1<:Real,N2<:Real,N3<:Real,N4<:Real,N5<:Real} <: AbstractMasterProblem{SimpleSubProblem} where {N1<:Number,N2<:Number,N3<:Number,N4<:Number,N5<:Number}
    costs::AbstractVector{N1}
    A::AbstractMatrix{N2}
    sense::AbstractVector{Char}
    b::AbstractVector{N3}
    vartypes::AbstractVector{Symbol}
    lb::AbstractVector{N4}
    ub::AbstractVector{N5}
    solver::MathProgBase.AbstractMathProgSolver
end

# TODO build convenient constructor functions for SimpleMasterProblem

end
