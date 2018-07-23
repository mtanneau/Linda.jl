"""
    Column

    Data structure for holding columns and related information.
"""
mutable struct Column 
    cost::Real          # native cost of the column
    col::AbstractVector # column coefficients
    isvertex::Bool      # whether column is vertex (true) or ray (false)
    isactive::Bool      # whether column is currently in the RMP
    problemidx::Int     # Index of sub-problem
end