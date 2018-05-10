"""
    SimpleProblem contains a full battery-included implementation of Linda Master and SubProblems
"""
module SimpleProblem

import MathProgBase

import Linda:
    AbstractSubProblem, AbstractMasterProblem, price,
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
function price(sp::SimpleSubProblem,π,σ, farkas_pricing = false)
    reduced_costs = [t[1] - t[2] for t in zip(sp.costs,π)] - σ
    result = MathProgBase.mixintprog(reduced_costs, sp.A, sp.sense, sp.b, sp.vartypes, sp.lb, sp.ub, sp.solver)
    final_status = find_status(result.status)
    return (final_status, [result.objval], hcat(result.sol))
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
