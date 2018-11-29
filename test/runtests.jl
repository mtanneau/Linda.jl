using Random
using LinearAlgebra
using Test


import MathProgBase
const MPB = MathProgBase

import Linda
import GLPKMathProgInterface:
    GLPKSolverLP, GLPKSolverMIP
import GLPK


const testdir = dirname(@__FILE__)

const test_files = [
    # include test file name here (without .jl extension)
    "status",
    "env",
    "master",
    "oracle_mip",
    "colgen"
]

for f in test_files
    tp = joinpath(testdir, "$(f).jl")
    include(tp)
end