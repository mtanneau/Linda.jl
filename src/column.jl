"""
    Column

    Data structure for holding columns and related information.
"""
struct Column
    cost::Float64  # Native cost

    spind::Int  # Index of convexity constraint

    # Column coefficients, in sparse form.
    # These do not include the convexity constraint
    rind::Vector{Int}
    rval::Vector{Float64}

    is_vertex::Bool # Whether column is vertex (true) or ray (false)
end