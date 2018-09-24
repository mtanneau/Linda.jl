using Random
using LinearAlgebra
using Test

import Linda
import GLPKMathProgInterface:
    GLPKSolverLP, GLPKSolverMIP
import MathProgBase
const MPB = MathProgBase

const testdir = dirname(@__FILE__)

const test_files = [
    # include test file name here (without .jl extension)
    "env",
    "master",
    "status",
    "oracle_mip",
    "colgen"
]

for f in test_files
    tp = joinpath(testdir, "$(f).jl")
    include(tp)
end