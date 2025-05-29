# __precompile__()

"""
Multi-Dimensional Bisection Method (MDBM) is an efficient and robust root-finding algorithm,
which can be used to determine whole high-dimensional sub-manifolds (points, curves, surfaces…)
of the roots of implicit non-linear equation systems, even in cases, where the number of unknowns
surpasses the number of equations.

# Examples
```julia
include("MDBM.jl")
using Reexport
@reexport using .MDBM

using PyPlot;
pygui(true);

function foo(x,y)
    x^2.0+y^2.0-2.0^2.0
end
function c(x,y) #only the c>0 domain is analysed
    x-y
end

ax1=Axis([-5,-2.5,0,2.5,5],"x") # initial grid in x direction
ax2=Axis(-5:2:5.0,"y") # initial grid in y direction

mymdbm=MDBM_Problem(foo,[ax1,ax2],constraint=c)
iteration=5 #number of refinements (resolution doubling)
solve!(mymdbm,iteration)


#points where the function foo was evaluated
x_eval,y_eval=getevaluatedpoints(mymdbm)

#interpolated points of the solution (approximately where foo(x,y) == 0 and c(x,y)>0)
x_sol,y_sol=getinterpolatedsolution(mymdbm)

fig = figure(1);clf()
scatter(x_eval,y_eval,s=5)
scatter(x_sol,y_sol,s=5)
```
"""
module MDBM
using StaticArrays
using LinearAlgebra

#using PrecompileTools: @setup_workload, @compile_workload    # this is a small dependency

export MDBM_Problem, Axis,
    solve!, interpolate!, refine!, checkneighbour!,
    axesextend!, getinterpolatedsolution, getevaluatedpoints, getevaluatedfunctionvalues, getevaluatedconstraintvalues,
    connect, triangulation, getinterpolatedgradient,
    interpsubcubesolution!, extract_paths

include("MDBM_types.jl")
include("MDBM_functions.jl")



#@setup_workload begin
#    # Putting some things in `@setup_workload` instead of `@compile_workload` can reduce the size of the
#    # precompile file and potentially make loading faster.
#    function foo_par3_codim1(x, y, z)
#        x^2.0 + y^2.0 + z^2.0 - 2.0^2.0
#    end
#    function c(x, y, z)
#        x^2.0 + y^2.0 - 0.5^2.0
#    end
#    @compile_workload begin
#        
#        mymdbm = MDBM_Problem(foo_par3_codim1, [-3.0:1.0, -1.0:3.0, -1.0:3.0], constraint=c)
#        solve!(mymdbm, 2, doThreadprecomp=false, verbosity=0)
#        xyz_sol = getinterpolatedsolution(mymdbm)
#    end
#end


end
