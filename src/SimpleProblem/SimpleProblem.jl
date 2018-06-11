"""
    SimpleProblem contains a full battery-included implementation of Linda Master and SubProblems
"""
module SimpleProblem

import MathProgBase
const MPB = MathProgBase

import Linda:
    AbstractSubProblem, AbstractMasterProblem, Column, PricingResult, MasterSolution,
    find_status, StatusError, StatusOptimal, StatusUnbounded, StatusInfeasible, ok, isinfeasible,
    compute_dual_variables!, get_subproblems, add_columns!, solve!,
    solve_pricing, getprobindex
    
export SimpleSubProblem, SimpleMasterProblem


include("./simpleSubProblem.jl")
include("./simpleMasterProblem.jl")

end  # module
