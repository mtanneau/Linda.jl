using Linda
using Cbc: CbcSolver
import Linda.SimpleProblem

@static if VERSION < v"0.7.0-DEV.2005"
    using Base.Test
else
    using Test
end

const testdir = dirname(@__FILE__)

const test_files = [
    "simple_subproblem",
    "cutting-stock",
    # include test file name here (without .jl extension)
]

for f in test_files
    tp = joinpath(testdir, "$(f).jl")
    include(tp)
end
