"""
    LindaOracleMIP

Solve a MIP, heuristically or to optimality, to generate a new column.
"""
mutable struct LindaOracleMIP <: AbstractLindaOracle
    index::Int  # Index of the sub-problem

    # Sub-problem data
    costs::AbstractVector{Float64}
    A_link::AbstractMatrix{Float64}
    A_sub::AbstractMatrix{Float64}
    row_lb::AbstractVector{Float64}
    row_ub::AbstractVector{Float64}
    vartypes::AbstractVector{Symbol}
    col_lb::AbstractVector{Float64}
    col_ub::AbstractVector{Float64}
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
        row_lb::AbstractVector{N4},
        row_ub::AbstractVector{N5},
        vartypes::AbstractVector{Symbol},
        col_lb::AbstractVector{N6},
        col_ub::AbstractVector{N7},
        solver::MPB.AbstractMathProgSolver
    ) where{N1<:Real,N2<:Real,N3<:Real,N4<:Real,N5<:Real,N6<:Real, N7<:Real}
        oracle = new()

        oracle.index = index

        oracle.costs = costs
        oracle.A_link = A_link
        oracle.A_sub = A_sub
        oracle.row_lb = row_lb
        oracle.row_ub = row_ub
        oracle.vartypes = vartypes
        oracle.col_lb = col_lb
        oracle.col_ub = col_ub
        oracle.solver = solver

        oracle.new_columns = Vector{Column}()
        oracle.status = Unknown
        oracle.pricing_dual_bound = -Inf

        return oracle
    end

end

function call_oracle!(
    oracle::LindaOracleMIP,
    π::AbstractVector{Tv},
    σ::Real;
    farkas=false
) where{Tv<:Real}

    # Dimension checks
    size(π, 1) == size(oracle.A_link, 1) || throw(DimensionMismatch("π has wrong size"))

    # Compute objective
    if farkas
        obj = -oracle.A_link' * π
    else
        obj = oracle.costs - oracle.A_link' * π
    end

    # Instanciate and solve sub-problem
    sp = MPB.LinearQuadraticModel(oracle.solver)
    MPB.loadproblem!(
        sp,
        oracle.A_sub,
        oracle.col_lb,
        oracle.col_ub,
        obj,
        oracle.row_lb,
        oracle.row_ub,
        :Min
    )
    MPB.setvartype!(sp, oracle.vartypes)
    MPB.optimize!(sp)

    status = findStatus(MPB.status(sp))
    oracle.status = status

    if status == Optimal
        # Sub-problem solved to optimality
        oracle.pricing_dual_bound = MPB.getobjbound(sp)
        x = MPB.getsolution(sp)
        oracle.new_columns = [
            Column(
                dot(oracle.costs, x),
                oracle.A_link * x,
                true,
                oracle.index
            )
        ]    
        
    elseif status == PrimalInfeasible
        # Sub-problem is infeasible
        oracle.new_columns = Vector{Column}()
        oracle.pricing_dual_bound = Inf

    elseif status == PrimalFeasible || status == PrimalDualFeasible
        # A primal solution is available
        oracle.pricing_dual_bound = MPB.getobjbound(sp)
        x = MPB.getsolution(sp)
        oracle.new_columns = [
            Column(
                dot(oracle.costs, x),
                oracle.A_link * x,
                true,
                oracle.index
            )
        ]

    elseif status == PrimalUnbounded

        #=
            TODO [Mathieu, 27/07/2018]
            
            Oracle should recover an unbounded ray as a new column with negative
            reduced cost.
            Some solvers (e.g. Cbc) do not implement that functionality, others
            may detect unboundedness but do not give access to an unbounded ray
            (e.g. if unboundedness is detected at presolve)
        =#
        warn("Unbounded sub-problem. Currently not supported.")
        oracle.status = Unknown
        
        #=
        # Sub-problem is unbounded
        # Return extreme ray
        ray = MPB.getunboundedray(sp)
        oracle.new_columns = [
            Column(
                dot(oracle.costs, ray),
                oracle.A_link * ray,
                false,
                oracle.index
            )
        ]
        oracle.pricing_dual_bound = -Inf
        =#

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

