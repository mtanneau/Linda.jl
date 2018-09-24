"""
    AbstractLindaParam{T}

Abstract representation of a parameter with valuetype `T`.
"""
abstract type AbstractLindaParam{T} end

"""
    get_param_name(::AbstractLindaParam)

Return name of parameter.
"""
function get_param_name end

"""
    get_param_type(::AbstractLindaParam)

Get parameter's value type.
"""
function get_param_type end

"""
    get_param_value(::AbstractLindaParam)

Return current value of parameter.
"""
function get_param_value end

"""
    set_param_default!(::AbstractLindaParam{T})

Set parameter to its default value.
"""
function set_param_default! end

"""
    set_param_value!(::AbstractLindaParam{T}, v::T)

Check if `v` is an admissible value for parameter and, if so, change parameter's
    value to `v`.
"""
function set_param_value! end

"""
    test_param_value(::AbstractLindaParam{T}, v::T)

Check whether value `v` is admissible for given parameter.
"""
function test_param_value end

"""
    LindaRealParam{T<:Real}

Container for numerical (real-valued) parameters.
"""
mutable struct LindaRealParam{T<:Real} <: AbstractLindaParam{T}
    name::Symbol  # Name of the parameter

    val::T  # Current parameter value
    min_val::T  # Minimum parameter value
    max_val::T  # Maximum parameter value
    def_val::T  # Default value

    LindaRealParam(name::Symbol, vdef::T, vmin::T, vmax::T) where{T<:Real} =
        new{T}(name, vdef, vmin, vmax, vdef)
    
end

const LindaIntParam = LindaRealParam{Int}
const LindaFloatParam = LindaRealParam{Float64}

get_param_name(par::LindaRealParam) = par.name

get_param_type(par::LindaRealParam{T}) where T = T

get_param_value(par::LindaRealParam{T}) where T = par.val

set_param_default!(par::LindaRealParam) = (par.val = par.def_val)

"""
    set_param_value!(p, v)

Set value of parameter `p` to `v`. Raises an error if `v` is not an admissible value.
"""
function set_param_value!(p::LindaRealParam, v::T) where{T<:Real}

    if test_param_value(p, v)
        p.val = v
    else
        error("$(p.name) must be between $(p.min_val) and $(p.max_val).")
    end

    return nothing
end

"""
    test_param_value(p, v)

Return whether `v` is an admissible value for parameter `p` 
"""
test_param_value(p::LindaRealParam, v::T) where{T<:Real} = (p.min_val <= v <= p.max_val)