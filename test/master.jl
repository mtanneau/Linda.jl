# Create initial RMP
m = 1  # number of linking constraints
R = 1  # number of sub-problems

b = ones(m)

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

# Generate some columns
for r in 1:R
    col = Linda.Column(1.0, ones(m), true, false, r, 0)
    Linda.add_columns!(mp, Set([col]))
    # Add column again, this should not do anything
    Linda.add_columns!(mp, Set([col]))
end

@test mp.num_columns_rmp == R
@test MPB.numvar(mp.rmp) == 2*m + R


# Solve RMP to optimality
Linda.solve_rmp!(mp)

@test mp.rmp_status == Linda.Optimal
@test mp.mp_status == Linda.PrimalFeasible
@test MPB.getobjval(mp.rmp) ≈ 1.0

# Infeasible RMP
# (de-activate artificial variables and )
MPB.setvarUB!(rmp, vcat(zeros(2*m), Inf*ones(R)))  # de-activate artificial vars
MPB.setconstrLB!(rmp, vcat(ones(R), -ones(m)))
MPB.setconstrUB!(rmp, vcat(ones(R), -ones(m)))
mp.rhs_constr_link = -ones(m)

Linda.solve_rmp!(mp)
@test mp.rmp_status == Linda.PrimalInfeasible  # RMP is primal infeasible
