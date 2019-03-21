

#include("MDBM.jl")
#using Reexport
#@reexport
using MDBM

using PyPlot;
pygui(true);


#-------------------


function foo(x,y)#,z)
    -4.0^2.0+x^2.0+y^2.0#+z^2
end
function c(x,y)#x,z)
    x-y#+z
end

ax1=Axis([-5,-2.5,0,2.5,5],"x")
ax2=Axis(-5:2:5.0,"b")
ax3=Axis(-5:2:5.0,"b")

mymdbm=MDBM_Problem(foo,[ax1,ax2])#ax3
iteration=2 #number of refinements (resolution doubling)
solve!(mymdbm,iteration,interpolationorder=1)

#evaluated points
P_eval=getevaluatedpoints(mymdbm)

#solution points
P_sol=getinterpolatedsolution(mymdbm)

fig = figure(1);clf()
scatter(P_eval[1],P_eval[2],s=2);
scatter(P_sol[1],P_sol[2],s=4);
#scatter3D(P_eval[1],P_eval[2],P_eval[3],s=2);
#scatter3D(P_sol[1],P_sol[2],P_sol[3],s=4);


myDT1=connect(mymdbm)
#plot solution lines one-by-one
fig = figure(11);clf()
for i in 1:length(myDT1)
    dt=myDT1[i]
    Ps=getinterpolatedsolution(mymdbm.ncubes[[dt...]],mymdbm)
    plot(Ps[1],Ps[2],marker=".")
    # plot3D(Ps[1],Ps[2],Ps[3],marker=".")
end





 using Profile


@time MDBM_Problem(foo,[ax1,ax2],constraint=c)
 Juno.@profiler  MDBM_Problem(foo,[ax1,ax2],constraint=c)

 mdbm_temp=MDBM_Problem(foo,[-5:5,-5:5]);
 @time solve!(mdbm_temp,5,interpolationorder=0)
 mdbm_temp=MDBM_Problem(foo,[-5:5,-5:5]);
 @time solve!(mdbm_temp,5,interpolationorder=1)

  mdbm_temp=MDBM_Problem(foo,[-5:5,-5:5]);
 Juno.@profiler  solve!(mdbm_temp,5,interpolationorder=0)
  mdbm_temp=MDBM_Problem(foo,[-5:5,-5:5]);
 Juno.@profiler  solve!(mdbm_temp,5,interpolationorder=1)
