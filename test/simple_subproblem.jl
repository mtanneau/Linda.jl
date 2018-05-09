import Linda.SimpleProblem

A = [
    1 1
    1 2
]

b = [3;3]
c = [-1;-1]
sense = ['<','<']
vartypes = [:Int,:Int]

sp = Linda.SimpleProblem.SimpleSubProblem(c,A,['<','<'],b, vartypes, [0,0],[Inf,Inf],CbcSolver())
(status, costs, columns) = Linda.solve(sp,[0.0;0.0],0.0)
@test status == Linda.StatusOptimal()
@test size(costs) == (1,)
@test size(columns) == (2,1)
@test columns[1,1] ≈ 3
@test columns[2,1] ≈ 0
