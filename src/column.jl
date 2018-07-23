"""
    Column

    Data structure for holding columns and related information.
"""
mutable struct Column 
    cost::Float64       # Native cost of the column
    col::AbstractVector{Float64} # column coefficients
    is_vertex::Bool      # Whether column is vertex (true) or ray (false)
    is_in_rmp::Bool      # Whether column is currently in the RMP
    idx_subproblem::Int  # Index of corresponding sub-problem
    idx_column::Int      # Index of the column in the RMP
end