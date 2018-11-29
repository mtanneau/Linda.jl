"""
    Status

These Status codes are based on those of MathOptInterface
- `PrimalFeasible`: A primal feasible solution has been found
- `DualFeasible` : A dual-feasible solution has been found
- `PrimalDualFeasible`: A primal- and a dual-feasible solutions have been found
- `Optimal`: A provably optimal solution has been found
- `PrimalInfeasible`: Primal is proved to be infeasible
- `PrimalUnbounded`: Problem is proved to be unbounded
- `Unknown`: No feasible solution nor proof of infeasibility yet
"""
@enum(Status,
    Unknown,
    PrimalFeasible,
    DualFeasible,
    PrimalDualFeasible,
    Optimal,
    PrimalInfeasible,
    PrimalUnbounded
)

Status(::T) where T = Unknown
Status(s::Symbol) = Status(Val(s))

Status(::Val{:Infeasible}) = PrimalInfeasible
Status(::Val{:Unbounded}) = PrimalUnbounded
Status(::Val{:Optimal}) = Optimal