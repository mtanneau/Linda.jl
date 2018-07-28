oracle = Linda.Oracle.LindaOracleMIP(
    1,
    [-1.0],
    ones(1, 1),
    ones(1, 1),
    [-Inf],
    [1.0],
    [:Cont],
    zeros(1),
    ones(1),
    ClpSolver()
)

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

mp = Linda.LindaMaster(R, m, b, ClpSolver(), oracle)
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
MPB.setvarUB!(mp.rmp, vcat(zeros(2*m), Inf*ones(mp.num_columns_rmp)))
# Change right-hand side
MPB.setconstrLB!(mp.rmp, vcat(ones(R), -ones(m)))
MPB.setconstrUB!(mp.rmp, vcat(ones(R), -ones(m)))
mp.rhs_constr_link = -ones(m)

Linda.solve_rmp!(mp)
@test mp.rmp_status == Linda.PrimalInfeasible  # RMP is primal infeasible