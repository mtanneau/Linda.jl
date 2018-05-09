"""
    AbstractStatus designs status signal returned after a problem is solved
    Basic ones are provided, other can be extended
"""
abstract type AbstractStatus end

"""
    StatusOptimal is returned if a problem could be solved to optimality
"""
struct StatusOptimal <: AbstractStatus end

"""
    StatusInfeasible is returned if a problem could not be solved
"""
struct StatusInfeasible <: AbstractStatus end

"""
    StatusUnbounded is returned if the solution in unbounded
"""
struct StatusUnbounded <: AbstractStatus end

"""
    StatusError is returned if something went wrong
    Use mesage and meta-information to document it
"""
struct StatusError <:AbstractStatus
    message::AbstractString
    meta::Dict{String,Any}
end

StatusError() = StatusError("", Dict{String,Any}())

find_status(::Type{Val{:Optimal}}) = StatusOptimal()
find_status(::Type{Val{:Infeasible}}) = StatusInfeasible()
find_status(::Type{Val{:Unbounded}}) = StatusUnbounded()

find_status(s::Symbol) = find_status(Val{s})

"""
    Construct status error from MathProgBase Symbol
"""
find_status(s) = StatusError("Other status from MathProgBase", Dict("status" => s))
