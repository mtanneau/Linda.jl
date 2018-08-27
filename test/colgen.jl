function add_initial_columns!(mp, m, R)

    # Add an initial set of columns
    Linda.add_columns!(
        mp,
        [Linda.Column(0.0, zeros(m), true, r) for r in 1:R]
    )

    return nothing
end

# Create initial RMP
srand(0)
n = 20  # original dimension
m = 1  # number of linking constraints
R = 10  # number of sub-problems
b = [R]  # Right-hand side of linking constraints

# Create oracle
# Sub-problem is a continuous knapsack
oracles = [
    Linda.Oracle.LindaOracleMIP(
        r,
        -collect(1:n),
        ones(1, n),
        ones(1, n),
        [-Inf],
        [10.0],
        [:Bin for _ in 1:n],
        zeros(n),
        ones(n),
        CbcSolver()
    )
    for r in 1:R
]


pool = Linda.Oracle.LindaOraclePool(oracles)

# Instanciate MP
mp = Linda.LindaMaster(
    R, m, b,
    ClpSolver()
)
add_initial_columns!(mp, m, R)

# Environment
env = Linda.LindaEnv()
env[:verbose] = 1  # Activate iteration log
env[:num_cgiter_max] = 10

# Solve problem with column generation
Linda.solve_colgen!(env, mp, pool)

# Make Master Infeasible, solve again
b_ = [n*R+1.0]
MPB.setconstrLB!(mp.rmp, vcat(ones(R), b_))
MPB.setconstrUB!(mp.rmp, vcat(ones(R), b_))
mp.rhs_constr_link = copy(b_)
MPB.setvarUB!(mp.rmp, vcat(zeros(2*m), Inf*ones(mp.num_columns_rmp)))
mp.dual_bound = -Inf

Linda.solve_colgen!(env, mp, pool)