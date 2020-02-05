using LinearAlgebra
using Test

import MathOptInterface
const MOI = MathOptInterface

import Linda
import GLPK

const testdir = dirname(@__FILE__)

@testset "Linda" begin
    include(joinpath(testdir, "env.jl"))
    include(joinpath(testdir, "master.jl"))
    include(joinpath(testdir, "colgen.jl"))
end