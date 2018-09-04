import Base: setindex!, getindex

include("params.jl")

"""
    LindaEnv

Optimizer environment. Holds parameters' values.
"""
mutable struct LindaEnv

    num_cols_rmp_max::LindaIntParam  # Maximum number of columns in RMP

    num_columns_max::LindaIntParam  # Maximum number of columns before pricing stops

    num_cgiter_max::LindaIntParam  # Maximum number of Column-Generation iterations

    time_limit::LindaFloatParam  # Time limit, in seconds

    tol_reduced_cost::LindaFloatParam  # Numerical tolerance for reduced costs

    verbose::LindaIntParam  # Verbosity level

    function LindaEnv()
        env = new()

        env.num_cols_rmp_max = LindaRealParam(:num_cols_rmp_max,
            typemax(Int64), zero(Int64), typemax(Int64)
        )

        env.num_columns_max = LindaRealParam(:num_columns_max,
            typemax(Int64), zero(Int64), typemax(Int64)
        )

        env.num_cgiter_max = LindaRealParam(:num_cgiter_max,
            typemax(Int64), zero(Int64), typemax(Int64)
        )

        env.time_limit = LindaRealParam(:time_limit, Inf, 0.0, Inf)

        env.tol_reduced_cost = LindaRealParam(:tol_reduced_cost, 10.0^-6, 0.0, Inf)

        env.verbose = LindaRealParam(:verbose, 0, 0, 1 )

        return env
    end

end

# TODO: use `getproperty` and `setproperty!` instead (> Julia 0.7)
getindex(env::LindaEnv, name::Symbol) = get_param_value(Core.getfield(env, name))

getindex(env::LindaEnv, ::Type{Val{:num_cols_rmp_max}}) = env.num_cols_rmp_max.val
getindex(env::LindaEnv, ::Type{Val{:num_columns_max}}) = env.num_columns_max.val
getindex(env::LindaEnv, ::Type{Val{:num_cgiter_max}}) = env.num_cgiter_max.val
getindex(env::LindaEnv, ::Type{Val{:time_limit}}) = env.time_limit.val
getindex(env::LindaEnv, ::Type{Val{:tol_reduced_cost}}) = env.tol_reduced_cost.val
getindex(env::LindaEnv, ::Type{Val{:verbose}}) = env.verbose.val

setindex!(env::LindaEnv, v, name::Symbol) = set_param_value!(Core.getfield(env, name), v)

"""
    setparam!(env; kwargs...)

Set parameters to specified values.
"""
function setparam!(env; kwargs...)
    for (s, v) in kwargs
        set_param_value!(Core.getfield(env, s), v)
    end
end

"""
    reset!(env)

Reset all parameters to their default values.
"""
function reset!(env)
    for s in fieldnames(env)
        set_param_default!(Core.getfield(env, s))
    end
    return nothing
end