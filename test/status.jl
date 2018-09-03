@test Linda.findStatus(:Infeasible) == Linda.PrimalInfeasible
@test Linda.findStatus(:Unbounded) == Linda.PrimalUnbounded
@test Linda.findStatus(:Optimal) == Linda.Optimal
@test Linda.findStatus(:Unknown) == Linda.Unknown

for i in 0:(length(instances(Linda.ProblemStatus))-1)
    @test Linda.findStatus(i) == Linda.ProblemStatus(i)
end