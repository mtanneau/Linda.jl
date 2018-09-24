# Default env constructor
env = Linda.LindaEnv()

# Test param functions
@test Linda.get_param_name(env.verbose) == :verbose
@test Linda.get_param_type(env.verbose) == Int64
@test Linda.get_param_type(env.time_limit) == Float64
@test Linda.get_param_value(env.verbose) == env.verbose.def_val
@test !Linda.test_param_value(env.verbose, -1)

# Test that `getindex(env, ::Type{Val{:###}}` have properly been set
for s in fieldnames(LindaEnv)
    @test env[Val{s}] == env[s]
end

# Check parameter assignment
env[:verbose] = 1.0
@test env.verbose.val == 1

env[:time_limit] = Int(1)
env[:time_limit] = 1.0
env[:time_limit] = BigFloat(1.0)
@test env.time_limit.val == 1.0

Linda.setparam!(env, verbose=1, time_limit=10)
@test env.verbose.val == 1
@test env.time_limit.val == 10.0

try 
    env[:time_limit] = -1.0
catch e
    @test e.msg == "time_limit must be between 0.0 and Inf."
end

Linda.reset!(env)
for s in fieldnames(LindaEnv)
    p = Core.getfield(env, s)
    @test env[s] == p.def_val
end