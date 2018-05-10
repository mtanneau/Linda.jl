"""
    SimpleProblem contains a full battery-included implementation of Linda Master and SubProblems
"""
module SimpleProblem

import MathProgBase

import Linda:
    AbstractSubProblem, AbstractMasterProblem, solve_pricing,
    find_status, StatusError, StatusOptimal, StatusUnbounded, StatusInfeasible

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
    price implementation for SimpleSubProblem
"""
function solve_pricing(sp::SimpleSubProblem, π, σ, farkas_pricing=false)

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
    final_status = find_status(result.status)

    # get solution
    if !ok(final_status)
        # Error when solving sub-problem: return no column
        return (final_status, Array{Float64,1}(), Matrix{Float64}(0, 0))
        # TODO: handle unbounded sub-problem
        #   In this case, a (primal) extreme ray is added to the master problem
    else
        col = result.sol  # column
        cost = dot(sp.costs, col)  # native cost
        rc = result.objval - σ  # reduced cost
    end

    # check that reduced cost is indeed negative
    if rc < -10.0^-6
        return (final_status, [cost], hcat(col))
    else
        # reduced cost is 0: return no column
        return (final_status, Array{Float64,1}(), Matrix{Float64}(0, 0))
    end

    # end of function
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
