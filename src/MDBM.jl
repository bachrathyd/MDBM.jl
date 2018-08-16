__precompile__()

module MDBM

export mdbm_problem, mdbm_object, refine!, checkncube!, checkneighbour!, interpolate!, DTconnect!
#export find_roots, refine_solution!

include("MDBM_types.jl")
include("MDBM_functions.jl")

end