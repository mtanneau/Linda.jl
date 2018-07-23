using Linda
import Cbc: CbcSolver
import Clp: ClpSolver
import MathProgBase
const MPB = MathProgBase

@static if VERSION < v"0.7.0-DEV.2005"
    using Base.Test
else
    using Test
end

const testdir = dirname(@__FILE__)

const test_files = [
    "simple_subproblem",
    "simple_masterproblem",
    # include test file name here (without .jl extension)
]

for f in test_files
    tp = joinpath(testdir, "$(f).jl")
    include(tp)
end
