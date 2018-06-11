"""
    SimpleMasterProblem

# Attributes
- `b::AbstractVector{T<:Real}`: right-hand side of linking constraints
- `rmp::MPB.AbstractLinearQuadraticModel`: Restricte Master Problem
- `sp::AbstractSubProblem`: Sub-Problem
"""
mutable struct SimpleMasterProblem{ST<:AbstractSubProblem, N1<:Real, V<:AbstractVector{N1}}<:AbstractMasterProblem{ST}
    
    #==================================================
        Master Problem
    ==================================================#
    numcon_link::Int    # Number of linking constraints
    b::V                # Right-hand side of linking constraints
    rmp::MPB.AbstractLinearQuadraticModel  # restricted Master Problem


    #==================================================
        Sub-Problem
    ==================================================#
    num_sp::Int     # Number of sub-problems
    sp::ST  # Sub-problems


    #==================================================
        Constructor
    ==================================================#
    SimpleMasterProblem(
        b::V,
        rmp::MPB.AbstractLinearQuadraticModel,
        sp::ST
    ) where {ST<:AbstractSubProblem, N1<:Real, V<:AbstractVector{N1}} =
        new{ST,N1,V}(size(b, 1), b, rmp, num_subproblems(sp), sp)
end

# TODO build convenient constructor functions for SimpleMasterProblem
function SimpleMasterProblem(
    b::V,
    senses::Vector{Char},
    solver::MPB.AbstractMathProgSolver,
    sp::ST
) where {N1<:Real, V<:AbstractVector{N1},ST<:AbstractSubProblem}

    # dimension check
    numcon_link = size(b, 1)  # number of linking constraints
    numcon_link == size(senses, 1) || DimensionMismatch("")
    num_sp = size(sp, 1)

    # instanciate rmp
    rmp = MPB.LinearQuadraticModel(solver)

    # create empty constraints with right-hand sides
    for j=1:numcon_link
        if senses[j] == '='
            MPB.addconstr!(rmp, [], [], b[j], b[j])
        elseif senses[j] == '<'
            MPB.addconstr!(rmp, [], [], -Inf, b[j])
        elseif senses[j] == '>'
            MPB.addconstr!(rmp, [], [], b[j], Inf)
        else
            error("Invalid input: senses[$(j)]=$(senses[j])")
        end
    end
    for j in 1:num_sp
        MPB.addconstr!(rmp, zeros(), zeros(), 1.0, 1.0)  # convexity constraint
    end
    
    return SimpleMasterProblem(numcon_link, b, rmp, num_sp, sp)
end


"""
    subproblem(::SimpleMasterProblem)

returns the SubProblem object
"""
function subproblem(mp::SimpleMasterProblem{ST}) where {ST<:AbstractSubProblem}
    return mp.sp
end

"""
    compute_dual_variables!(::SimpleMasterProblem)
    
Solve the RMP to optimality and return new dual iterate.
"""
function compute_dual_variables!(mp::SimpleMasterProblem{ST}) where {ST<:AbstractSubProblem}

    # solve RMP
    MPB.optimize!(mp.rmp)
    rmp_status = find_status(MPB.status(mp.rmp))

    # check RMP status
    nlinkingconstrs = size(mp.b, 1)
    if !ok(rmp_status)
        # unexpected status when solving RMP
        # TODO: handle dual unboundedness in RMP
        π = zeros(0,)
        σ = zeros(1,)
    else
        # RMP solved to optimality
        dualsolution = MPB.getconstrduals(mp.rmp)
        π = dualsolution[1:nlinkingconstrs]
        σ = dualsolution[nlinkingconstrs+1:end]
    end

    return MasterSolution(rmp_status, π, σ)
end

function add_columns!(mp::SimpleMasterProblem{ST}, columns::Vector{Column}) where ST<:AbstractSubProblem
    ncols = size(columns, 1)
    ncolsadded = 0

    constridx = collect(1:mp.numcon_link)

    for column in columns
        # add column (assumed to have negative reduced cost)
        constridx_ = [constridx ; mp.numcon_link+column.problemidx]
        
        if column.isactive
            continue
        end
        ncolsadded += 1
        if column.isvertex
            # extreme vertex
            constrcoeff = vcat(column.col, [1.0])
        else
            # extreme ray
            constrcoeff = vcat(column.col, [0.0])
        end
        MPB.addvar!(mp.rmp, constridx_, constrcoeff, 0.0, Inf, column.cost)
        column.isactive = true
    end
    return ncolsadded
end