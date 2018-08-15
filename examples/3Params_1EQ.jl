include((x->x[1:find(t->t=='/', x)[end]])(pwd())*"src/MDBM.jl")
using MDBM
using Plots
plotly()
# plotlyjs()


# -----3 parameter , 1 implicit equation ---------------
f1(x,y,z) = [x*x+y*y-z*z-1]
ax=[collect(-5.0:5.0),collect(-5.0:5.0),-5:5]

mdbm_solutionSurf=mdbm_problem(f1,ax,fullrefinenum=2)# run all the steps automatcally

scatter(mdbm_solutionSurf.posinterp[1,:],mdbm_solutionSurf.posinterp[2,:],mdbm_solutionSurf.posinterp[3,:], leg=false, xlim=(-5,5), ylim=(-5,5), zlim=(-5,5), markersize=1)
