"""
    AbstractColumn

    Abstract structure for holding columns.
"""

mutable struct AbstractColumn end

"""
    getnativecost
    Return the (native) cost of a column
"""
function getnativecost(::AbstractColumn)
    warn("Implement getcolcost for concrete implementation!")
    return 0.0
end

