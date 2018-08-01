import Base: setindex!, getindex

include("linda_params.jl")

"""
    LindaEnv

Optimizer environment. Holds parameters' values.
"""
mutable struct LindaEnv
    params::Dict{Symbol, AbstractLindaParam}  # Dictionnary of parameters

    #===========================================================================
            Default parameters
    ===========================================================================#
    function LindaEnv()
        env = new()
        env.params = Dict{Symbol, AbstractLindaParam}()
        
        # Maximum number of columns in RMP
        add_param!(env, LindaNumParam(
            :num_cols_rmp_max,
            typemax(Int64),
            zero(Int64),
            typemax(Int64)
        ))

        # Maximum number of Column-Generation iterations
        add_param!(env, LindaNumParam(
            :num_cgiter_max,
            typemax(Int64),
            zero(Int64),
            typemax(Int64)
        ))

        # Time-limit
        add_param!(env, LindaNumParam(
            :time_limit,
            Inf,
            0.0,
            Inf
        ))

        # Reduced-cost tolerance
        # Should be no smaller than teh RMP solver's reduced cost tolerance,
        #   otherwise cycling will occur.
        add_param!(env, LindaNumParam(
            :tol_reduced_cost,
            10.0^-6,
            0.0,
            Inf
        ))

        # Verbosity level
        add_param!(env, LindaNumParam(
            :verbose,
            0,
            0,
            1
        ))
        return env
    end



end


function getindex(env::LindaEnv, name::Symbol)
    
    if !haskey(env.params, name)
        error("Parameter $name does not exist.")
    else
        return get_param_value(env.params[name])
    end
end

function setindex!(env::LindaEnv, v, name::Symbol)

    if !haskey(env.params, name)
        error("Parameter $name does not exist.")
    else
        return set_param_value!(env.params[name], v)
    end
end

"""
    add_param!(env::LindaEnv, param::AbstractLindaParam)

Add parameter to the environment. If there exists a parameter with same name,
    it will be replaced.
"""
function add_param!(env::LindaEnv, param::AbstractLindaParam)
    name = get_param_name(param)

    # Add parameter / Replace if already exists.
    env.params[name] = param
end

