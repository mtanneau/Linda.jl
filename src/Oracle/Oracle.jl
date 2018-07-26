module Oracle

import MathProgBase
const MPB = MathProgBase

import Linda:
    Column,
    ProblemStatus,
    Unknown, PrimalFeasible, DualFeasible, PrimalDualFeasible, Optimal, PrimalInfeasible, PrimalUnbounded,
    findStatus


"""
    AbstractLindaOracle

Each call to the oracle will generate a finite number of columns.

An oracle must return:
- The set of newly generated columns, if any
- The best known lower bound, if applicable
"""
abstract type AbstractLindaOracle end


"""
    call_oracle!

Execute a call to the oracle.
"""
function call_oracle! end


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


include("oracle_mip.jl")  # Sub-problem is solved as a MIP

end  # Oracle module