# __precompile__()

module MDBM
using Reexport
@reexport using StaticArrays
@reexport using LinearAlgebra

export MDBM_Problem, Axis,
 solve!, interpolate!, refine!, checkneighbour!,
 axesextend!, getinterpolatedsolution, getevaluatedpoints, getevaluatedfunctionvalues, getevaluatedconstraintvalues,
 connect, triangulation

include("MDBM_types.jl")
include("MDBM_functions.jl")

end
