"""
    Column

    Data structure for holding columns and related information.
"""
mutable struct Column{Tv<:Real, Tc<:AbstractVector{Tv}}
    cost::Float64        # Native cost of the column
    col::Tc              # column coefficients
    is_vertex::Bool      # Whether column is vertex (true) or ray (false)
    is_in_rmp::Bool      # Whether column is currently in the RMP
    idx_subproblem::Int  # Index of corresponding sub-problem
    idx_column::Int      # Index of the column in the RMP

    Column(
        cost,
        col::Tc,
        is_vertex,
        idx_subproblem
    ) where{Tv<:Real, Tc<:AbstractArray{Tv}} = new{Tv, Tc}(cost, col, is_vertex, false, idx_subproblem, -1)
end

"""
    get_reduced_cost()

Compute reduced cost of column given current dual iterate
"""
get_reduced_cost(
    c::Column{Tv, Tc},
    π::AbstractVector{Float64},
    σ::Real,
    farkas=false
) where{Tv<:Real, Tc<:AbstractArray{Tv}} = (!farkas)*c.cost - dot(π, c.col) - σ

get_reduced_cost(
    c::Column{Tv, Tc},
    π::AbstractVector{Float64},
    σ::AbstractVector{Float64},
    farkas=false
) where{Tv<:Real, Tc<:AbstractArray{Tv}} = 
    get_reduced_cost(c, π, σ[c.idx_subproblem], farkas)