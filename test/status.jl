@test Linda.Status(:Infeasible) == Linda.PrimalInfeasible
@test Linda.Status(:Unbounded) == Linda.PrimalUnbounded
@test Linda.Status(:Optimal) == Linda.Optimal
@test Linda.Status(:Unknown) == Linda.Unknown