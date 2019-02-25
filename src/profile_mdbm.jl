

include("MDBM.jl")
using Reexport
@reexport using .MDBM

using PyPlot;
pygui(true);


#-------------------






function foo(x,y)
    x^2.0+y^2.0-4.0^2.0
end
function c(x,y)
    x-y
end

ax1=Axis([-5,-2.5,0,2.5,5],"x")
ax2=Axis(-5:2:5.0,"b")

mymdbm=MDBM_Problem(foo,[ax1,ax2],constraint=c)
iteration=4 #number of refinements (resolution doubling)
solve!(mymdbm,iteration)

#evaluated points
x_eval,y_eval=getevaluatedpoints(mymdbm)

#solution points
x_sol,y_sol=getinterpolatedsolution(mymdbm)

fig = figure(1);clf()
scatter(x_eval,y_eval,s=2)
scatter(x_sol,y_sol,s=4);


using Profile



Juno.@profiler  mymdbm=MDBM_Problem(foo,[ax1,ax2],constraint=c)
Juno.@profiler solve!(mymdbm,iteration)
