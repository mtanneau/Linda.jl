n = 3
m = 1
k = 2.0
A_link = ones(1, n)
b = [k]
c1 = [-1.0, -2.0, -3.0]
c2 = [-1.5, -2.5, -3.5]

sp1 = Linda.SimpleProblem.SimpleSubProblem(
    c1,
    A_link,
    Matrix{Float64}(0, n),
    Vector{Char}(0, ),
    Vector{Float64}(0,),
    [:Bin for i=1:n],
    zeros(n), ones(n),
    CbcSolver(),
    1
)

sp2 = Linda.SimpleProblem.SimpleSubProblem(
    c2,
    A_link,
    Matrix{Float64}(0, n),
    Vector{Char}(0, ),
    Vector{Float64}(0,),
    [:Bin for i=1:n],
    zeros(n), ones(n),
    CbcSolver(),
    2
)

subproblems = [sp1, sp2]

# create MasterProblem
# instanciate RMP
rmp = MathProgBase.LinearQuadraticModel(ClpSolver())

# add linking GUB constraint
MathProgBase.addconstr!(rmp, [], [], -Inf, k)  # empty constraints (no variable)

# add convexity constraint
for sp in subproblems
    MathProgBase.addconstr!(rmp, [], [], 1.0, 1.0)

    # initialize master with zero-valued column
    # in our case this column is feasible
    MathProgBase.addvar!(rmp, [collect(1:m); (m+sp.problemidx)], [0., 1.], 0., Inf, 0.)
end

mp = Linda.SimpleProblem.SimpleMasterProblem(b, rmp, subproblems);

Linda.solve!(mp)

obj_ = MPB.getobjval(rmp)
@test abs(obj_ - (-6.5)) <= 10.0^-8
