__precompile__()

module MDBM

export mdbm_problem, mdbm_object, refine!, checkncube!, checkneighbour!, interpolate!, DTconnect!
#export find_roots, refine_solution!

include("MDBM_types.jl")
include("MDBM_functions.jl")

end

function f(x::Float64,y::Float64,z::Float64)::Float64
    x^2+y^2+z^2-4.0
end
ccc=typeof(f(1.0,2.0,34.0))
typeof(f)

a=5
aff =MemoizationFunction(f)
typeof(f1)
aff =MemoizationFunction(f,)

returntype(sin,(Int64,))
listss=sort([(2,3),(3,3),(2,5)])
Lia=in((3,3),listss)
findfirst(x -> x==(3,3),((2,3),(3,3),(2,5)))

fun=MemoizationFunction(x->x^3)
fun.f(5)

bb=Any[]
fun(12)
