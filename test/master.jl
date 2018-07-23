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
    col = Linda.Column(0.0, ones(m), true, false, r, 0)
    Linda.add_columns!(mp, Set([col]))
end

Linda.solve_rmp!(mp)