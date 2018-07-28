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
R = 2  # number of sub-problems
b = [5.0]  # Right-hand side of linking constraints

# Create oracle
# Sub-problem is a continuous knapsack
oracles = [Linda.Oracle.LindaOracleMIP(
    r,
    -collect(1:n),
    ones(1, n),
    ones(1, n),
    [-Inf],
    [10.0],
    [:Cont for _ in 1:n],
    zeros(n),
    ones(n),
    CbcSolver()
) for r in 1:R]

handler = Linda.Oracle.LindaOracleHandler(oracles)

mp = Linda.LindaMaster(R, m, b, ClpSolver(), handler)
add_initial_columns!(mp, m, R)

Linda.solve_colgen!(mp)