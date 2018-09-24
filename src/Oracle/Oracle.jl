module Oracle

import MathProgBase
const MPB = MathProgBase

import Linda:
    Column, get_reduced_cost,
    Status, Unknown, PrimalFeasible, DualFeasible, PrimalDualFeasible, Optimal, PrimalInfeasible, PrimalUnbounded

export
    AbstractLindaOracle, call_oracle!, get_new_columns, get_sp_dual_bound


"""
    AbstractLindaOracle

Each call to the oracle will generate a finite number of columns.

An oracle must return:
- The set of newly generated columns, if any
- The best known lower bound, if applicable
"""
abstract type AbstractLindaOracle end


"""
    query!

Execute a query to the oracle.
"""
function query! end


"""
    get_oracle_status

Get status of oracle.
"""
function get_oracle_status end


"""
    get_num_new_columns

Get number of newly generated columns.
"""
function get_num_new_columns end


"""
    get_new_columns

Retrieve newly generated columns.
"""
function get_new_columns end


"""
    get_sp_dual_bound

Get best known dual bound from last call to the oracle. If no dual bound is 
    available (e.g. sub-problem solved heuristically), infinity is returned.
"""
function get_sp_dual_bound end

"""
    get_sp_best_objective

Get best objective, i.e. most negative reduced cost
"""
function get_sp_best_objective end

include("oracle_pool.jl")
include("oracle_mip.jl")  # Solve sub-problem as a MIP

end  # Oracle module