# Create a knapsack sub-problem for oracle tests

n = 1  # number of items
c = [5.0]  # values of items
w = 5.0*ones(1, 1)  # weights of items, as a 1-row matrix
v = 10  # total capacity
solver_mip = CbcSolver()

# Bounded sub-problem
# Min -x
# s.t. x <= 1
#      x binary
oracle_bounded = Linda.Oracle.LindaOracleMIP(
    1,
    [-1.0],
    ones(1, 1),
    ones(1, 1),
    [-Inf],
    [1.0],
    [:Bin],
    zeros(1),
    ones(1),
    solver_mip
)

# Unbounded sub-problem
# Min -x
# s.t. x >= 1
#      x integer
oracle_unbound = Linda.Oracle.LindaOracleMIP(
    1,
    [-1.0],
    ones(1, 1),
    ones(1, 1),
    [1.0],
    [Inf],
    [:Cont],
    zeros(1),
    Inf * ones(1),
    solver_mip
)

# Min -x
# s.t. x = 2
#      x binary
oracle_infeas = Linda.Oracle.LindaOracleMIP(
    1,
    [-1.0],
    ones(1, 1),
    ones(1, 1),
    [2.0],
    [2.0],
    [:Bin],
    zeros(1),
    ones(1),
    solver_mip
)

# Dual variables
π = [10.0]
σ = 0.0

# Bounded oracle
Linda.Oracle.call_oracle!(oracle_bounded, π, σ)  # solve sub-problem
status = Linda.Oracle.get_oracle_status(oracle_bounded)
if status == Linda.Optimal
    # at least one column should have been generated
    @test Linda.Oracle.get_num_new_columns(oracle_bounded) >= 1
    cols = Linda.Oracle.get_new_columns(oracle_bounded)
    col = cols[1]
    # check value of reduced costs and dual bound
    @test Linda.get_reduced_cost(col, π, σ) ≈ -1.0 - π[1]
    @test Linda.Oracle.get_sp_dual_bound(oracle_bounded) ≈ -1.0 - π[1]
    # check computation of column's cost
    @test col.cost ≈ -1.0

elseif status == Linda.PrimalFeasible
    # at least one column should have been generated
    @test Linda.Oracle.get_num_new_columns(oracle_bounded) >= 1
    cols = Linda.Oracle.get_new_columns(oracle_bounded)
    col = cols[1]
    # check value of reduced costs and dual bound
    @test Linda.get_reduced_cost(col, π, σ) >= -1.0 - π[1]
    # check computation of column's cost
    @test Linda.Oracle.get_sp_dual_bound(oracle_bounded) <= -1.0 - π[1]

elseif status == Linda.Unknown
    warn("Bounded oracle: Solver should have solved problem")
else
    error("")
end

# Unbounded sub-problem
Linda.Oracle.call_oracle!(oracle_unbound, π, σ)
status = Linda.Oracle.get_oracle_status(oracle_unbound)
if status == Linda.PrimalUnbounded
    @test Linda.Oracle.get_num_new_columns(oracle_unbound) >= 1
    cols = Linda.Oracle.get_new_columns(oracle_unbound)
    col = cols[1]
    @test Linda.get_reduced_cost(col, π, σ) ≈ -1.0 - π[1]
    @test Linda.Oracle.get_sp_dual_bound(oracle_unbound) ≈ -Inf
    @test col.cost ≈ -1.0

elseif status == Linda.PrimalFeasible
    @test Linda.Oracle.get_num_new_columns(oracle_unbound) >= 1
    cols = Linda.Oracle.get_new_columns(oracle_unbound)
    col = cols[1]
    @test Linda.get_reduced_cost(col, π, σ) >= -1.0 - π[1]
    @test Linda.Oracle.get_sp_dual_bound(oracle_unbound) <= -1.0 - π[1]

elseif status == Linda.Unknown
    warn("Unbounded returned Unknown status")
else
    error("")
end

# Infeasible sub-problem
Linda.Oracle.call_oracle!(oracle_infeas, π, σ)
status = Linda.Oracle.get_oracle_status(oracle_infeas)
@test typeof(status) == Linda.ProblemStatus