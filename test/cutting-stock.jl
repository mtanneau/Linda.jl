const prices = [14, 65, 95, 90, 75,160,150, 95,190,120,145,200,175,200,210,250]
const widths = [10, 30, 45, 48, 50, 60, 69, 75, 78, 84, 90, 96,100,120,144,150]
const demand = [22, 50,156, 88,108,157, 43, 68,120,180, 72, 19,162, 33, 70,  5]
const rollcost = 750
const maxwidth = 225
const NW = length(widths)

initial_cols = [i!=j?0: floor(Int,maxwidth//widths[i]) for i=1:NW,j=1:NW]

knapsack_problem = SimpleProblem.SimpleSubProblem(
    -zeros(NW,), hcat(widths)', ['<'], [maxwidth],
    [:Int for _ in 1:NW],
    zeros(Int,NW,),
    [floor(Int,maxwidth/widths[i]) for i in 1:NW],
    CbcSolver()
)

mutable struct CuttingWidthMaster{SP<:SimpleProblem.SimpleSubProblem} <: Linda.AbstractMasterProblem{SP}
    current_cols::Matrix{Int64}
    current_costs::Vector{Float64}
    demand::Vector{Int64}
    sp::SP
    solver::ClpSolver
end

# implement interface
Linda.subproblem(cw::CuttingWidthMaster) = cw.sp

function Linda.add_columns!(cw::CuttingWidthMaster,columns::Vector{Linda.Column})
    for col in columns
        cw.current_cols = hcat(cw.current_cols,[round(Int,a) for a in col.col])
        push!(cw.current_costs,1.0)
    end
    return length(columns)
end

function Linda.compute_dual_variables!(cw::CuttingWidthMaster)
    nwidths = size(cw.current_cols)[1]
    result = MathProgBase.linprog(
        cw.current_costs, cw.current_cols, '>', cw.demand, cw.solver
    )
    status = Linda.find_status(result.status)
    return Linda.MasterSolution(status, -result.attrs[:lambda], [1.0])
end

cw = CuttingWidthMaster(initial_cols,ones(NW,),demand,knapsack_problem,ClpSolver())

include("ref_columns.jl")
result_status = Linda.solve!(cw, maxcols = 1000)
@test result_status == Linda.StatusOptimal()
@test all(excess_rolls.==ref_excess)

# master_problem = SimpleProblem.SimpleMasterProblem(
#     initial_cols, ['>' for _ in 1:NW], demand, ClpSolver(), knapsack_problem
# )

# result_status = solve!(master_problem, maxcols = 1000)
# @test result_status == Linda.StatusOptimal()
# cols = master_problem.columns
# include("ref_columns.jl")
# @test all(cols.==patterns)

# excess_rolls = round(Int,master_problem.columns * master_problem.solution - demand)
# @test all(excess_rolls.==ref_excess)
