# Create initial RMP
rmp = GLPK.Optimizer()
m = 2  # number of linking constraints
R = 2  # number of sub-problems
b = zeros(m)  # Right-hand side of linking constraints

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
        MOI.EqualTo(b[i])
    )
    push!(con_link, cidx)
end

# Add linking variables
vlink = MOI.VariableIndex[]

# Add initial columns
cols = MOI.VariableIndex[]
for r in 1:R
    vidx, _ = MOI.add_constrained_variable(rmp, MOI.GreaterThan(0.0))
    push!(cols, vidx)

    # Set objective coefficient
    MOI.modify(rmp,
        MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(),
        MOI.ScalarCoefficientChange(vidx, 1.0)
    )
end


mp = Linda.Master(rmp, con_cvx, con_link, b, vlink, cols)

# TODO: Add columns