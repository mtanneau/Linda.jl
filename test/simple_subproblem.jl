@test Linda.find_status(:Optimal) == Linda.StatusOptimal()
@test Linda.find_status(:Infeasible) == Linda.StatusInfeasible()
@test Linda.find_status(:Unbounded) == Linda.StatusUnbounded()

A = [
    1 1
    1 2
]

b = [3;3]
c = [-1;-1]
sense = ['<','<']
vartypes = [:Int,:Int]

sp = Linda.SimpleProblem.SimpleSubProblem(c,A,['<','<'],b, vartypes, [0,0],[Inf,Inf],CbcSolver())
presult = Linda.solve_pricing(sp,[0.0;0.0],[0.0])
@test presult.status == Linda.StatusOptimal()
@test size(presult.columns) == (1,)
# @test size(columns) == (2,1)
# @test columns[1,1] ≈ 3
# @test columns[2,1] ≈ 0
