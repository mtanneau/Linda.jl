const prices = [14, 65, 95, 90, 75,160,150, 95,190,120,145,200,175,200,210,250]
const widths = [10, 30, 45, 48, 50, 60, 69, 75, 78, 84, 90, 96,100,120,144,150]
const demand = [22, 50,156, 88,108,157, 43, 68,120,180, 72, 19,162, 33, 70,  5]
const rollcost = 750
maxwidth = 225

knapsack_problem = SimpleProblem.SimpleSubProblem(
    -prices, hcat(widths)', ['<'], [maxwidth],
    [:Int for _ in 1:length(widths)], 
    spzeros(Int,length(prices)),
    [floor(Int,maxwidth/widths[i]) for i in 1:length(prices)],
    CbcSolver()
)

master_problem = SimpleProblem.SimpleMasterProblem{SimpleProblem.SimpleSubProblem}(
    A, senses, b, CbcSolver(), knapsack_problem
)

    A::AbstractMatrix{N1},
    senses::AbstractVector{Char},
    b::AbstractVector{N2},
    solver::MathProgBase.AbstractMathProgSolver,
    sp::ST

result_status = solve!(master_problem, maxcols = 1000)
@test result_status == Linda.StatusOptimal()
cols = master_problem.columns
include("ref_columns.jl")
@test all(cols.==patterns)

excess_rolls = round(Int,master_problem.columns * master_problem.solution - demand)
@test all(excess_rolls.==ref_excess)
