using LinearAlgebra

import GLPK

# Needed until gets merged in GLPK.jl
MOI.get(::GLPK.Optimizer, ::MOI.BarrierIterations) = 0
MOI.get(::GLPK.Optimizer, ::MOI.SimplexIterations) = 0

# Dummy oracle: minimize c'x over the hypercube
mutable struct DummyOracle <: Linda.Oracle
    id::Int
    dim::Int

    c::Vector{Float64}
    x::Vector{Float64}

    farkas::Bool
    π::Vector{Float64}
    σ::Float64

    DummyOracle(id::Int, dim::Int, c) = new(id, dim, c, zeros(dim), false, zeros(dim), 0.0)

end

function Linda.update!(o::DummyOracle, farkas, π, σ)
    o.farkas = farkas
    o.π .= π
    o.σ = σ
    return nothing
end

Linda.optimize!(::DummyOracle) = nothing

function Linda.get_columns(o::DummyOracle)
    o.x .= 0.0
    rc = -o.σ

    if o.farkas
        # Farkas pricing
        for i in 1:o.dim
            if -o.π[i] <= 0.0
                o.x[i] = 1.0
                rc -= o.π[i]
            end
        end
    else
        # Regular pricing
        for i in 1:o.dim
            if o.c[i] - o.π[i] < 0.0
                o.x[i] = 1.0
                rc += o.c[i] - o.π[i]
            end
        end
    end

    col = Linda.Column(dot(o.x, o.c), o.id, collect(1:o.dim), copy(o.x), true)
    return [(col, rc)]
end

Linda.get_dual_bound(o::DummyOracle) = o.farkas ? -Inf : (dot(o.x, o.c) - dot(o.x, o.π) - o.σ)

Linda.get_objective_value(o::DummyOracle) = o.farkas*dot(o.x, o.c) - dot(o.x, o.π) - o.σ

# Create initial RMP
import GLPK
rmp = GLPK.Optimizer()
# Disable RMP solver's output
MOI.set(rmp, MOI.Silent(), true)
# Set objective sense to minimization, otherwise all the duals are wrong
MOI.set(rmp, MOI.ObjectiveSense(), MOI.MIN_SENSE)


m = 4  # number of linking constraints
R = 1  # number of sub-problems
b = [0.2, 0.4, 0.6, 0.8]

# Create the convexity constraints
con_cvx = MOI.ConstraintIndex[]
for r in 1:R
    cidx = MOI.add_constraint(
        rmp,
        MOI.ScalarAffineFunction(MOI.ScalarAffineTerm{Float64}[], 0.0),
        MOI.EqualTo(1.0)
    )
    push!(con_cvx, cidx)
end
# Create linking constraints
con_link = MOI.ConstraintIndex[]
for i in 1:m
    cidx = MOI.add_constraint(
        rmp,
        MOI.ScalarAffineFunction(MOI.ScalarAffineTerm{Float64}[], 0.0),
        MOI.LessThan(b[i])
    )
    push!(con_link, cidx)
end

mp = Linda.Master(rmp, con_cvx, con_link, b, MOI.VariableIndex[], MOI.VariableIndex[])

# Add a single column
# Linda.add_column!(mp,
#     Linda.Column(
#         0.0, 1, Int[], Float64[], true
# ))

# Add a single column
# Linda.add_column!(mp,
#     Linda.Column(
#         -4.0, 1, [1], [1.0], true
# ))
# Linda.add_column!(mp,
#     Linda.Column(
#         -8.0, 1, [2], [1.0], true
# ))
# Linda.add_column!(mp,
#     Linda.Column(
#         -12.0, 1, [3], [1.0], true
# ))
# Linda.add_column!(mp,
#     Linda.Column(
#         -16.0, 1, [4], [1.0], true
# ))
# The following three give an optimal solution
# Linda.add_column!(mp,
#     Linda.Column(
#         -16.0, 1, [1, 3], [1.0, 1.0], true
# ))
# Linda.add_column!(mp,
#     Linda.Column(
#         -40.0, 1, [1, 2, 3, 4], [1.0, 1.0, 1.0, 1.0], true
# ))

# Setup the sub-problems
oracles = [DummyOracle(1, m, [-5.0, -4.0, -3.0, -2.0])]

# CG parameters
env = Linda.LindaEnv()
env.num_cgiter_max.val = 10
env.verbose.val = 1

# Solve the problem with Column Generation
Linda.solve_colgen!(env, mp, oracles)