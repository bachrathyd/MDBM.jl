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
# scatter(x_eval,y_eval,s=5)
# scatter(x_sol,y_sol,s=5)
```
"""
module MDBM
using StaticArrays
using LinearAlgebra
using FunctionWrappers

using PrecompileTools: @setup_workload, @compile_workload

export MDBM_Problem, Axis,
    solve!, interpolate!, refine!, checkneighbour!,
    axesextend!, getinterpolatedsolution, getevaluatedpoints, getevaluatedfunctionvalues, getevaluatedconstraintvalues,
    connect, triangulation, getinterpolatedgradient,
    interpsubcubesolution!, extract_paths,
    recreate
include("MDBM_types.jl")
include("MDBM_functions.jl")



@setup_workload begin
#     # using Plots - for testing the output
# 
#     # ------------ 2D ------------
#     function foo_par2_codim1(x, y)
#         x^2.0 + y^2.0 - 2.0^2.0
#     end
#     function coo_par2(x, y)
#         x + y
#     end
#     # ------------ 3D ------------
#     function foo_par3_codim1(x, y, z)
#         x^2.0 + y^2.0 + z^2.0 - 2.0^2.0
#     end
#     function foo_par3_codim2(x, y, z)
#         return x^2.0 + y^2.0 + z^2.0 - 2.0^2.0, x - sin(z * 2)
#     end
#     function coo_par3(x, y, z)
#         x + y + z
#     end
#     # ------------ 4D ------------
#     function foo_par4_codim2(x, y, z, w)
#         return x^2.0 + y^2.0 + z^2.0 + w^2.0 - 2.0^2.0, x - sin(z * 2)
#     end
#     function foo_par4_codim3(x, y, z, w)
#         return x^2.0 + y^2.0 + z^2.0 + w^2.0 - 2.0^2.0, x - sin(z * 2), y - cos(w * 2)
#     end
#     function coo_par4(x, y, z, w)
#         x + y + z + w
#     end
# 
# 
     @compile_workload begin
#         
#         # ------------ 2D ------------
#         sin(foo_par2_codim1(1.2,3.21))
#         mymdbm = MDBM_Problem(foo_par2_codim1, [-3.0:3.0, -3.0:3.0])
#         solve!(mymdbm, 1, doThreadprecomp=false, verbosity=1)
#         mymdbm = MDBM_Problem(foo_par2_codim1, [-3.0:3.0, -3.0:3.0])
#         solve!(mymdbm, 1, doThreadprecomp=true, verbosity=1)
#         xy_sol = getinterpolatedsolution(mymdbm)
#         # scatter(xy_sol...,)
#         mymdbm = MDBM_Problem(foo_par2_codim1, [-3.0:3.0, -3.0:3.0], constraint=coo_par2)
#         solve!(mymdbm, 1, doThreadprecomp=false, verbosity=1)
#         mymdbm = MDBM_Problem(foo_par2_codim1, [-3.0:3.0, -3.0:3.0], constraint=coo_par2)
#         solve!(mymdbm, 1, doThreadprecomp=true, verbosity=1)
#         xy_sol = getinterpolatedsolution(mymdbm)
#         DT1 = connect(mymdbm)
#         # scatter(xy_sol...,)
# 
#         # ------------ 3D ------------
#         mymdbm = MDBM_Problem(foo_par3_codim1, [-3.0:3.0, -3.0:3.0, -3.0:3.0])
#         solve!(mymdbm, 1, doThreadprecomp=false, verbosity=1)
#         mymdbm = MDBM_Problem(foo_par3_codim1, [-3.0:3.0, -3.0:3.0, -3.0:3.0])
#         solve!(mymdbm, 1, doThreadprecomp=true, verbosity=1)
#         xyz_sol = getinterpolatedsolution(mymdbm)
#         # scatter(xyz_sol...,)
#         mymdbm = MDBM_Problem(foo_par3_codim1, [-3.0:3.0, -3.0:3.0, -3.0:3.0], constraint=coo_par3)
#         solve!(mymdbm, 1, doThreadprecomp=false, verbosity=1)
#         mymdbm = MDBM_Problem(foo_par3_codim1, [-3.0:3.0, -3.0:3.0, -3.0:3.0], constraint=coo_par3)
#         solve!(mymdbm, 1, doThreadprecomp=true, verbosity=1)
#         xyz_sol = getinterpolatedsolution(mymdbm)
#         DT1 = connect(mymdbm)
#         # scatter(xyz_sol...,)
# 
#         mymdbm = MDBM_Problem(foo_par3_codim2, [-3.0:3.0, -3.0:3.0, -3.0:3.0])
#         solve!(mymdbm, 1, doThreadprecomp=false, verbosity=1)
#         mymdbm = MDBM_Problem(foo_par3_codim2, [-3.0:3.0, -3.0:3.0, -3.0:3.0])
#         solve!(mymdbm, 1, doThreadprecomp=true, verbosity=1)
#         xyz_sol = getinterpolatedsolution(mymdbm)
#         # scatter(xyz_sol...,)
#         mymdbm = MDBM_Problem(foo_par3_codim2, [-3.0:3.0, -3.0:3.0, -3.0:3.0], constraint=coo_par3)
#         solve!(mymdbm, 1, doThreadprecomp=false, verbosity=1)
#         mymdbm = MDBM_Problem(foo_par3_codim2, [-3.0:3.0, -3.0:3.0, -3.0:3.0], constraint=coo_par3)
#         solve!(mymdbm, 1, doThreadprecomp=true, verbosity=1)
#         xyz_sol = getinterpolatedsolution(mymdbm)
#         DT1 = connect(mymdbm)
#         # scatter(xyz_sol...,)
# 
#         # ------------ 4D ------------
#         mymdbm = MDBM_Problem(foo_par4_codim2, [-3.0:3.0, -3.0:3.0, -3.0:3.0, -3.0:3.0])
#         solve!(mymdbm, 1, doThreadprecomp=false, verbosity=1)
#         mymdbm = MDBM_Problem(foo_par4_codim2, [-3.0:3.0, -3.0:3.0, -3.0:3.0, -3.0:3.0])
#         solve!(mymdbm, 1, doThreadprecomp=true, verbosity=1)
#         xyzw_sol = getinterpolatedsolution(mymdbm)
#         # scatter(xyzw_sol[1:3]...,)
#         mymdbm = MDBM_Problem(foo_par4_codim2, [-3.0:3.0, -3.0:3.0, -3.0:3.0, -3.0:3.0], constraint=coo_par4)
#         solve!(mymdbm, 1, doThreadprecomp=false, verbosity=1)
#         mymdbm = MDBM_Problem(foo_par4_codim2, [-3.0:3.0, -3.0:3.0, -3.0:3.0, -3.0:3.0], constraint=coo_par4)
#         solve!(mymdbm, 1, doThreadprecomp=true, verbosity=1)
#         xyzw_sol = getinterpolatedsolution(mymdbm)
#         DT1 = connect(mymdbm)
#         # scatter(xyzw_sol[1:3]...,)
# 
#         mymdbm = MDBM_Problem(foo_par4_codim3, [-3.0:3.0, -3.0:3.0, -3.0:3.0, -3.0:3.0])
#         solve!(mymdbm, 1, doThreadprecomp=false, verbosity=1)
#         mymdbm = MDBM_Problem(foo_par4_codim3, [-3.0:3.0, -3.0:3.0, -3.0:3.0, -3.0:3.0])
#         solve!(mymdbm, 1, doThreadprecomp=true, verbosity=1)
#         xyzw_sol = getinterpolatedsolution(mymdbm)
#         DT1 = connect(mymdbm)
#         # scatter(xyzw_sol[1:3]...,)
# 
#         mymdbm = MDBM_Problem(foo_par4_codim3, [-3.0:3.0, -3.0:3.0, -3.0:3.0, -3.0:3.0], constraint=coo_par4)
#         solve!(mymdbm, 1, doThreadprecomp=false, verbosity=1)
#         mymdbm = MDBM_Problem(foo_par4_codim3, [-3.0:3.0, -3.0:3.0, -3.0:3.0, -3.0:3.0], constraint=coo_par4)
#         solve!(mymdbm, 1, doThreadprecomp=true, verbosity=1)
#         xyzw_sol = getinterpolatedsolution(mymdbm)
#         DT1 = connect(mymdbm)
#         # scatter(xyzw_sol[1:3]...,)
 include("../deps/precompile_workload.jl")
   end
end


end
