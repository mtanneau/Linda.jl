const prices = [14, 65, 95, 90, 75,160,150, 95,190,120,145,200,175,200,210,250]
const widths = [10, 30, 45, 48, 50, 60, 69, 75, 78, 84, 90, 96,100,120,144,150]
const demand = [22, 50,156, 88,108,157, 43, 68,120,180, 72, 19,162, 33, 70,  5]
const rollcost = 750
const maxwidth = 225
const NW = length(widths)

initial_cols = [i!=j?0: floor(Int,maxwidth//widths[i]) for i=1:NW,j=1:NW]

knapsack_problem = SimpleProblem.SimpleSubProblem(
    -prices, hcat(widths)', ['<'], [maxwidth],
    [:Int for _ in 1:NW], 
    spzeros(Int,NW,),
    [floor(Int,maxwidth/widths[i]) for i in 1:NW],
    CbcSolver()
)

master_problem = SimpleProblem.SimpleMasterProblem{SimpleProblem.SimpleSubProblem}(
    initial_cols, ['>' for _ in 1:NW], b, CbcSolver(), knapsack_problem
)

result_status = solve!(master_problem, maxcols = 1000)
@test result_status == Linda.StatusOptimal()
cols = master_problem.columns
include("ref_columns.jl")
@test all(cols.==patterns)

excess_rolls = round(Int,master_problem.columns * master_problem.solution - demand)
@test all(excess_rolls.==ref_excess)
