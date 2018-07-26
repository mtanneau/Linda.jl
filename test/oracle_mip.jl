# Create a knapsack sub-problem for oracle tests

n = 10  # number of items
c = rand(1:10, n)  # values of items
w = rand(1:20, 1, n)  # weights of items, as a 1-row matrix
v = 10*n  # total capacity

solver_mip = CbcSolver()

oracle = Linda.Oracle.LindaOracleMIP(
    1,
    -c,
    ones(1, n),
    w,
    ['<'],
    [v],
    [:Bin for _ in 1:n],
    zeros(n),
    ones(n),
    solver_mip
)

π = [0.0]
σ = 0.0

Linda.Oracle.call_oracle!(oracle, π, σ)

status = Linda.Oracle.get_oracle_status(oracle)
@test typeof(status) == Linda.ProblemStatus

if status == Linda.Optimal
    @test Linda.Oracle.get_num_new_columns(oracle) >= 1
    cols = Linda.Oracle.get_new_columns(oracle)
    col = cols[1]
    @test Linda.get_reduced_cost(col, π, σ) == Linda.Oracle.get_sp_dual_bound(oracle)

elseif status == Linda.PrimalUnbounded
    @test Linda.Oracle.get_num_new_columns(oracle) >= 1
end