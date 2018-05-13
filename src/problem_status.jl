"""
    AbstractStatus designs status signal returned after a problem is solved.
    Basic ones are provided, other can be defined by inheritance and defining
        the `ok()` function.
"""
abstract type AbstractStatus end

"""
    ok is a function used to check if a status is fine
        (i.e., iterations of the column generation should continue)
"""
ok(::AbstractStatus) = true

"""
    StatusOptimal is returned if a problem is solved to proven optimality.
"""
struct StatusOptimal <: AbstractStatus end

ok(::StatusOptimal) = true

"""
    StatusInfeasible is returned if a problem is proven to be infeasible.
    In the case of Linear Programs, the solver should also provide
        a dual extreme ray as proof of infeasibility.
"""
struct StatusInfeasible <: AbstractStatus end
ok(::StatusInfeasible) = false

"""
    StatusFeasible is returned when a feasible solution is found, but optimality
        was not proven.
"""
struct StatusFeasible <: AbstractStatus end
ok(::StatusFeasible) = false  # to avoid cycling in the column-generation

"""
    StatusUnbounded is returned if the problem is proven to be unbounded.
    In the case of Linear Programs, the solver should also provide a primal
        extreme ray as proof of unboundedness.
"""
struct StatusUnbounded <: AbstractStatus end
ok(::StatusUnbounded) = false

"""
    StatusTimeout is returned if no solution was found within certain
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

find_status(s::Symbol)::AbstractStatus = find_status(Val{s})

"""
    Construct status error from MathProgBase Symbol
"""
find_status(s) = StatusError("Other status from MathProgBase", Dict("status" => s))
