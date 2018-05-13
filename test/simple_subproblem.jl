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

sp = SimpleProblem.SimpleSubProblem(c,A,['<','<'],b, vartypes, [0,0],[Inf,Inf],CbcSolver())
(status, costs, columns) = Linda.solve_pricing(sp,[0.0;0.0],[0.0])
@test status == Linda.StatusOptimal()
@test size(costs) == (1,)
@test size(columns) == (2,1)
@test columns[1,1] ≈ 3
@test columns[2,1] ≈ 0

# type stability tests
@inferred Linda.solve_pricing(sp,[0.0;0.0],0.0)
