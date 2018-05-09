import Linda.SimpleProblem

A = [
    1 1
    1 2
]

b = [3;3]
c = [1;1]
sense = ['<','<']
vartypes = [:Int,:Int]

sp = Linda.SimpleProblem.SimpleSubProblem(c,A,['<','<'],b, vartypes, [0,0],[Inf,Inf],CbcSolver())
(status, costs, columns) = Linda.solve(sp,[0.0;0.0],0.0)
@test status == Linda.StatusOptimal()
