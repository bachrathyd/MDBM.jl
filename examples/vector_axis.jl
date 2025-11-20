using MDBM
using PyPlot
pygui(true)

using LinearAlgebra
### Axis with vector valued nodes

# Task: find the solution along the eges of a triangle


# the function which defines the region on the plane of the triangle
fig = figure(1);clf()
Circ_mdbm=MDBM_Problem((x,y,r)->x*x+y*y-r*r,[-2.0:2.0,-2.0:2.0,0.75:0.2:1.25])
solve!(Circ_mdbm,3, doThreadprecomp=false, verbosity=0);
x_sol,y_sol=getinterpolatedsolution(Circ_mdbm)
scatter(x_sol,y_sol,s=2);

# Computing along the edge of a triangle:
# The ticks of the ax1 is the nodes of the triangle
ax1=Axis([[-1.0,1.0],[-1.0,-1.0],[1.5,-0.5],[-1.0,1.0]],"v")
ax2=Axis(0.75:0.2:1.25,"r")

#plotting the tirangle
for k in 1:length(ax1.ticks)-1
plot([ax1.ticks[k][1],ax1.ticks[k+1][1]],[ax1.ticks[k][2],ax1.ticks[k+1][2]])
end

function foo_Vect_Float(v,r)
   norm(v)-r
end

Vect_mdbm=MDBM_Problem(foo_Vect_Float,[ax1,ax2])
# the tree nodes of the triangle is not a proper initial mesh, so ax1 it is refined 3 times
refine!(Vect_mdbm,directions=[1,1,1])
scatter([P[1] for P in Vect_mdbm.axes[1]],[P[2] for P in Vect_mdbm.axes[1]],s=20)
solve!(Vect_mdbm,6)

#solution points
V_sol,r_sol=getinterpolatedsolution(Vect_mdbm) #Note V_sol is a Vector of Vectors

scatter([P[1] for P in V_sol],[P[2] for P in V_sol],s=6);

#Plotting the result in 3D along the first and second element of the vector (from ax1) + the redial coordinate (from ax2)
fig = figure(2);clf()
scatter3D([P[1] for P in V_sol],[P[2] for P in V_sol],r_sol,s=6);
