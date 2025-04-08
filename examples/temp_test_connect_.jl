using MDBM
using Plots
#using PyPlot
#pygui(true)



function foo(x,y,z)
    #x^2.0+y^2.0+1.0/sin(z)^2.0-4.0^2.0
    x^2.0+y^2.0+2*sin(z)^2.0-4.0^2.0,x-sin(z)
end
function c(x,y,z)
    -x+y+sin(z)
end

ax1=Axis([-5,-2.5,00.1,2.5,5],"x")
ax2=Axis(-5:2:5.1,"y")
ax3=Axis(-5:2.02:5.1,"z")

mymdbm=MDBM_Problem(foo,[ax1,ax2,ax3])#,constraint=c)
mymdbm=MDBM_Problem(foo,[ax1,ax2,ax3],constraint=c)
iteration=3 #number of refinements (resolution doubling)
solve!(mymdbm,iteration)

#evaluated points
x_eval,y_eval,z_eval=getevaluatedpoints(mymdbm)

#solution points
x_sol,y_sol,z_sol=getinterpolatedsolution(mymdbm)

#fig = figure(1);clf()
#scatter(x_eval,y_eval,z_eval,ms=1)
scatter(x_sol,y_sol,z_sol,ms=2)


myDT1=connect(mymdbm)



plot()
for i in 1:length(myDT1)
    dt=myDT1[i]
    P1x,P1y,P1z=getinterpolatedsolution(mymdbm.ncubes[dt[1]],mymdbm)
    P2x,P2y,P2z=getinterpolatedsolution(mymdbm.ncubes[dt[2]],mymdbm)
    plot!([P1x,P2x],[P1y,P2y],[P1z,P2z])
end 
plot!()


#DT2=triangulation(myDT1)


