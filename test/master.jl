function create_mp(m, R, b)

    rmp = MPB.LinearQuadraticModel(ClpSolver())
    # Convexity constraints
    for r in 1:R
        MPB.addconstr!(rmp, Vector{Float64}(), Vector{Float64}(), 1.0, 1.0)
    end

    # linking constraints
    for i in 1:m
        MPB.addconstr!(rmp, Vector{Float64}(), Vector{Float64}(), b[i], b[i])
    end

    # Artificial variables
    for i in 1:m
        MPB.addvar!(rmp, [R+i], [1.0], 0.0, Inf, 10^4)  # slack
        MPB.addvar!(rmp, [R+i], [-1.0], 0.0, Inf, 10^4)  # surplus
    end

    # Instanciate RMP
    mp = Linda.LindaMaster(rmp, R, m, b)

    return mp
end

function add_initial_columns!(mp, m, R)

    # Add an initial set of columns
    Linda.add_columns!(
        mp,
        [Linda.Column(0.0, zeros(m), true, r) for r in 1:R]
    )

    # Add some other columns
    for r in 1:R
        col = Linda.Column(1.0, rand(m), true, r)
        Linda.add_column!(mp, col)
        # Add column again, this should not do anything
        Linda.add_columns!(mp, [col])
    end

    return nothing
end

# Create initial RMP
srand(0)
m = 2  # number of linking constraints
R = 2  # number of sub-problems
b = zeros(m)  # Right-hand side of linking constraints

mp = create_mp(m, R, b)
add_initial_columns!(mp, m, R)

@test mp.num_columns_rmp == 2*R
@test MPB.numvar(mp.rmp) == 2*m + 2*R


# Solve RMP to optimality
Linda.solve_rmp!(mp)
s = MPB.getreducedcosts(mp.rmp)
# check computation of reduced costs
for i in 1:mp.num_columns_rmp
    @test s[2*m+i] ≈ Linda.get_reduced_cost(mp.active_columns[i], mp.π, mp.σ)
end

@test mp.rmp_status == Linda.Optimal
@test mp.mp_status == Linda.PrimalFeasible
@test MPB.getobjval(mp.rmp) ≈ 0.0

# Make RMP infeasible
# de-activate artificial variables
MPB.setvarUB!(rmp, vcat(zeros(2*m), Inf*ones(mp.num_columns_rmp)))
# Change right-hand side
MPB.setconstrLB!(rmp, vcat(ones(R), -ones(m)))
MPB.setconstrUB!(rmp, vcat(ones(R), -ones(m)))
mp.rhs_constr_link = -ones(m)

Linda.solve_rmp!(mp)
@test mp.rmp_status == Linda.PrimalInfeasible  # RMP is primal infeasible
=#