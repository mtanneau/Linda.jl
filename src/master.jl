"""
    Master

"""
mutable struct Master

    #===========================================================================
        RMP data
    ===========================================================================#
    
    rmp::MOI.ModelLike  # Restricted Master Problem

    # Constraints information
    con_cvx::Vector{MOI.ConstraintIndex}  # Convexity constraints in RMP
    con_link::Vector{MOI.ConstraintIndex}  # Linking constraints in RMP

    rhs_link::Vector{Float64}  # Right-hand side of linking constraints

    π::Vector{Float64}  # Shadow prices (linking constraints)
    σ::Vector{Float64}  # Shadow marginal costs (convexity constraints)

    # Columns information
    var_link::Vector{MOI.VariableIndex}  # Linking variables (may be artificial)
    columns::Vector{MOI.VariableIndex}  # Indices of RMP variables
    

    #===========================================================================
        Column-Generation
    ===========================================================================#

    #=
        Master Problem Status

        The Master Problem is:
        - Feasible if a feasible solution has been found in the RMP.
        - Optimal if solved to optimality, i.e. there exists a primal
            feasible solution, and pricing fails to identify a column of
            negative reduced cost.
        - Infeasible if proven infeasible, i.e. RMP is infeasible and
            Farkas pricing fails to cut the infeasibility ray.
        - Unbounded if the RMP is unbounded
    =#
    mp_status::MOI.TerminationStatusCode  # Status of MasterProblem
    mp_gap::Float64

    primal_lp_bound::Float64  # Primal bound for the linear (DW) relaxation
    primal_ip_bound::Float64  # Primal bound for the integer problem
    dual_bound::Float64   # Lagrange dual bound
    dual_bound_estimate::Float64  # Estimate of the Lagrange dual bound


    #===========================================================================
        Constructor
    ===========================================================================#

    function Master(
        rmp::MOI.ModelLike,
        con_cvx::Vector{MOI.ConstraintIndex},
        con_link::Vector{MOI.ConstraintIndex},
        rhs_link::Vector{Float64},
        var_link::Vector{MOI.VariableIndex},
        columns::Vector{MOI.VariableIndex}
    )
        # Check indices
        mcvx = length(con_cvx)
        mlink = length(con_link)
        for cidx in con_cvx
            MOI.is_valid(rmp, cidx) || error("Invalid constraint index $cidx")
        end
        for cidx in con_link
            MOI.is_valid(rmp, cidx) || error("Invalid constraint index $cidx")
        end
        length(rhs_link) == mlink || error(
            "There are $mlink constraints but RHS has length $(length(rhs_link))"
        )

        # Check variable indices
        for vidx in var_link
            MOI.is_valid(rmp, vidx) || error("Invalid variable index $vidx")
        end
        for vidx in columns
            MOI.is_valid(rmp, vidx) || error("Invalid variable index $vidx")
        end

        # Instanciate model
        mp = new(
            rmp,
            con_cvx, con_link, rhs_link,
            zeros(mlink), zeros(mcvx),
            var_link, columns,
            MOI.OPTIMIZE_NOT_CALLED, Inf, Inf, Inf, -Inf, -Inf
        )

        return mp
    end
end


"""
    add_column!(mp, col)

Add column `col` to Master Problem `mp`.

All columns are non-negative by default.
"""
function add_column!(mp::Master, col::Column)

    # add column to RMP
    vidx, _ = MOI.add_constrained_variable(mp.rmp, MOI.GreaterThan(0.0))

    # Set objective coefficient
    MOI.modify(mp.rmp,
        MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(),
        MOI.ScalarCoefficientChange(vidx, col.cost)
    )

    # Set column coefficient for convexity constraint
    if col.is_vertex
        MOI.modify(mp.rmp, mp.con_cvx[col.spind], MOI.ScalarCoefficientChange(vidx, 1.0))
    end
    # Set column coefficients for linking constraints
    for (i, val) in zip(col.rind, col.rval)
        MOI.modify(mp.rmp, mp.con_link[i], MOI.ScalarCoefficientChange(vidx, val))
    end

    # Book-keeping
    push!(mp.columns, vidx)

    # Done
    return vidx
end