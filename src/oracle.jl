"""
    Oracle

Each call to the oracle will generate a finite number of columns.

An oracle must return:
- The set of newly generated columns, if any
- The best known lower bound, if applicable
"""
abstract type Oracle end

"""
    update!(o::Oracle, π, σ)

Update objective function with new dual variables
"""
function update! end

"""
    optimize!(o::Oracle)

Solve the sub-problem.
"""
function optimize! end


"""
    get_columns

Return newly generated columns and corresponding reduced costs.
"""
function get_columns end


"""
    get_dual_bound

Get best known dual bound from last call to the oracle. If no dual bound is 
    available (e.g. sub-problem solved heuristically), infinity is returned.
"""
function get_dual_bound end

"""
    get_objective_value

Get best objective value of best solution, i.e., the most negative reduced cost.
"""
function get_objective_value end