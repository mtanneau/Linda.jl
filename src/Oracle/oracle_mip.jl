"""
    LindaOracleMIP

Solve a MIP, heuristically or to optimality, to generate a new column.
"""
mutable struct LindaOracleMIP{N1<:Real,N2<:Real,N3<:Real,N4<:Real,N5<:Real,N6<:Real} <: AbstractLindaOracle
    index::Int  # Index of the sub-problem

    # Sub-problem data
    costs::AbstractVector{N1}
    A_link::AbstractMatrix{N2}
    A_sub::AbstractMatrix{N3}
    senses::Vector{Char}
    b::AbstractVector{N4}
    vartypes::AbstractVector{Symbol}
    lb::AbstractVector{N5}
    ub::AbstractVector{N6}
    solver::MPB.AbstractMathProgSolver

    # Result of last oracle call
    # Includes newly generated columns and associated data
    new_columns::Vector{Column}
    status::ProblemStatus  # status of sub-problem
    pricing_dual_bound::Float64  # dual bound for sub-problem

    # TODO: Column pool (?)

    function LindaOracleMIP(
        index::Int,
        costs::AbstractVector{N1},
        A_link::AbstractMatrix{N2},
        A_sub::AbstractMatrix{N3},
        senses::AbstractVector{Char},
        b::AbstractVector{N4},
        vartypes::AbstractVector{Symbol},
        lb::AbstractVector{N5},
        ub::AbstractVector{N6},
        solver::MPB.AbstractMathProgSolver
    ) where{N1<:Real,N2<:Real,N3<:Real,N4<:Real,N5<:Real,N6<:Real}
        oracle = new{N1, N2, N3, N4, N5, N6}()

        oracle.index = index

        oracle.costs = costs
        oracle.A_link = A_link
        oracle.A_sub = A_sub
        oracle.senses = senses
        oracle.b = b
        oracle.vartypes = vartypes
        oracle.lb = lb
        oracle.ub = ub
        oracle.solver = solver

        oracle.new_columns = Vector{Column}()
        oracle.status = Unknown
        oracle.pricing_dual_bound = -Inf

        return oracle
    end

end

function call_oracle!(
    oracle::LindaOracleMIP{N1, N2, N3, N4, N5, N6},
    π::AbstractVector{Tv},
    σ::Real;
    farkas=false
) where{N1<:Real, N2<:Real, N3<:Real, N4<:Real, N5<:Real, N6<:Real, Tv<:Real}

    # Dimension checks
    size(π, 1) == size(oracle.A_link, 1) || throw(DimensionMismatch())

    # Compute objective
    if farkas
        obj = -oracle.A_link' * π
    else
        obj = oracle.costs - oracle.A_link' * π
    end

    # Solve sub-problem
    result = MPB.mixintprog(
        obj,
        oracle.A_sub,
        oracle.senses,
        oracle.b,
        oracle.vartypes,
        oracle.lb,
        oracle.ub,
        oracle.solver
    )

    status = findStatus(result.status)
    oracle.status = status

    if status == PrimalInfeasible
        oracle.new_columns = Vector{Column}()
        oracle.pricing_dual_bound = Inf

    elseif status == PrimalUnbounded
        oracle.new_columns = [
            Column(
                dot(oracle.costs, result.sol),
                oracle.A_link * result.sol,
                false,
                oracle.index
            )
        ]
    elseif status == Optimal
        oracle.pricing_dual_bound = result.attrs[:objbound]
        oracle.new_columns = [
            Column(
                dot(oracle.costs, result.sol),
                oracle.A_link * result.sol,
                true,
                oracle.index
            )
        ]
    elseif status == PrimalFeasible
        oracle.new_columns = [
            Column(
                dot(oracle.costs, result.sol),
                oracle.A_link * result.sol,
                true,
                oracle.index
            )
        ]
    else
        # TODO
        warn("Pricing status $(status) not handled")
    end

    return nothing
end

get_oracle_status(oracle::LindaOracleMIP) = oracle.status

get_num_new_columns(oracle::LindaOracleMIP) = size(oracle.new_columns, 1)

get_new_columns(oracle::LindaOracleMIP) = oracle.new_columns

get_sp_dual_bound(oracle::LindaOracleMIP) = oracle.pricing_dual_bound

