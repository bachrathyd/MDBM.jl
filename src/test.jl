using MDBM

using PyPlot;
pygui(true);

#------------




using StaticArrays
using LinearAlgebra
include("MDBM_types.jl")
include("MDBM_functions.jl")

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

#------------------------Svec test

a = SVector{7}([-3,-2,-1,0,1, 2, 3])
mymdbm=MDBM_Problem((x)->x^2-2,[a])
solve!(mymdbm,5)
println(getinterpolatedsolution(mymdbm))

#------------------------Svec test


a = Axis([-3,-2,-1,0,1, 2, 3])
b = Axis([-3,-2,-1,0,1, 2, 3])
axS = @SVector [a,b]
axV = [a,b]

axes1=[Axis(ax) for ax in axS]
axes2=[Axis(ax) for ax in axV]

mymdbm=MDBM_Problem((x,y)->x^2-y,axS)
solve!(mymdbm,2)
#println(getinterpolatedsolution(mymdbm))
size(mymdbm.ncubes)


mymdbm=MDBM_Problem((x,y)->x^2-y,axV)
solve!(mymdbm,2)
#println(getinterpolatedsolution(mymdbm))
size(mymdbm.ncubes)
