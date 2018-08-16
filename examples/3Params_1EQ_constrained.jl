include((x->x[1:find(t->t=='/', x)[end]])(pwd())*"src/MDBM.jl")
using MDBM
using Plots
plotly()
# plotlyjs()

# -----3 parameter , 1 implicit equation, 1 implicit equation for constrain ---------------
f1(x,y,z) = [x*x+y*y-z*z]
fc(x,y,z) = y+sin(x)+2
 
mdbm_solutionConstrainedSurface=mdbm_problem(f1,ax,fullrefinenum=2,fconstrain=fc)
scatter(mdbm_solutionConstrainedSurface.posinterp[1,:],mdbm_solutionConstrainedSurface.posinterp[2,:],mdbm_solutionConstrainedSurface.posinterp[3,:], leg=false, xlim=(-5,5), ylim=(-5,5), zlim=(-5,5), markersize=1)

