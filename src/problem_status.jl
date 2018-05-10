"""
    AbstractStatus designs status signal returned after a problem is solved
    Basic ones are provided, other can be defined by inheritance and defining ok()
"""
abstract type AbstractStatus end

"""
    ok is a function used to check if a status is
    fine (iterations of the column generation should continue)
"""
ok(::AbstractStatus) = true

"""
    StatusOptimal is returned if a problem could be solved to optimality
"""
struct StatusOptimal <: AbstractStatus end

ok(::StatusOptimal) = true

"""
    StatusInfeasible is returned if a problem could not be solved
"""
struct StatusInfeasible <: AbstractStatus end
ok(::StatusInfeasible) = false

"""
    StatusUnbounded is returned if the solution in unbounded
"""
struct StatusUnbounded <: AbstractStatus end
ok(::StatusUnbounded) = false

"""
    StatusTimeout is returned if the solution was not found under certain
    computation limits (number of columns, time, etc)
"""
struct StatusTimeout <: AbstractStatus end
ok(::StatusTimeout) = false

"""
    StatusError is returned if something went wrong
    Use mesage and meta-information to document it
"""
struct StatusError <:AbstractStatus
    message::AbstractString
    meta::Dict{String,Any}
end
ok(::StatusError) = false

StatusError() = StatusError("", Dict{String,Any}())

find_status(::Type{Val{:Optimal}}) = StatusOptimal()
find_status(::Type{Val{:Infeasible}}) = StatusInfeasible()
find_status(::Type{Val{:Unbounded}}) = StatusUnbounded()

find_status(s::Symbol) = find_status(Val{s})

"""
    Construct status error from MathProgBase Symbol
"""
find_status(s) = StatusError("Other status from MathProgBase", Dict("status" => s))
