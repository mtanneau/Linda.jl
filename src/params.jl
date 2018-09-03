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
    LindaNumParam{T<:Real}

Container for numerical (real-valued) parameters.
"""
mutable struct LindaNumParam{T<:Real} <: AbstractLindaParam{T}
    name::Symbol  # Name of the parameter

    val::T  # Current parameter value
    min_val::T  # Minimum parameter value
    max_val::T  # Maximum parameter value
    def_val::T  # Default value

    LindaNumParam(name::Symbol, vdef::T, vmin::T, vmax::T) where{T<:Real} =
        new{T}(name, vdef, vmin, vmax, vdef)
    
end

const LindaIntParam = LindaNumParam{Int}
const LindaFloatParam = LindaNumParam{Float64}

get_param_name(par::LindaNumParam) = par.name

get_param_type(par::LindaNumParam{T}) where T = T

get_param_value(par::LindaNumParam{T}) where T = par.val

set_param_default!(par::LindaNumParam) = (par.val = par.def_val)

function set_param_value!(par::LindaNumParam, v::T) where{T<:Real}

    if !test_param_value(par, v)
        error("$(par.name) must be between $(par.min_val) and $(par.max_val).")
    else
        par.val = v
    end

    return nothing
end

test_param_value(par::LindaNumParam, v::T) where{T<:Real} = (par.min_val <= v <= par.max_val)