"""
    ProblemStatus

These status codes are based on those of MathOptInterface
- `PrimalFeasible`: A primal feasible solution has been found
- `DualFeasible` : A dual-feasible solution has been found
- `PrimalDualFeasible`: A primal- and a dual-feasible solutions have been found
- `Optimal`: A provably optimal solution has been found
- `PrimalInfeasible`: Primal is proved to be infeasible
- `PrimalUnbounded`: Problem is proved to be unbounded
- `Unknown`: No feasible solution nor proof of infeasibility yet
"""
@enum ProblemStatus Unknown PrimalFeasible DualFeasible PrimalDualFeasible Optimal PrimalInfeasible PrimalUnbounded

findStatus(::Symbol) = Unknown
findStatus(::Type{Val{:Infeasible}}) = PrimalInfeasible
findStatus(::Type{Val{:Unbounded}}) = PrimalUnbounded
findStatus(::Type{Val{:Optimal}}) = Optimal