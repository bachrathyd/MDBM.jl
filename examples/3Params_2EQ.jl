include((x->x[1:find(t->t=='/', x)[end]])(pwd())*"src/MDBM.jl")
using MDBM
using Plots
plotly()
# plotlyjs()

# -----3 parameter , 2 implicit equations ---------------
f2(x,y,z) = [x*x+y*y-z*z-1,y+sin(x)+2]
ax=[collect(-5.0:5.0),collect(-5.0:5.0),-5:5]

@time mdbm_solutionLine=mdbm_problem(f2,ax,fullrefinenum=2)

scatter!(mdbm_solutionLine.posinterp[1,:],mdbm_solutionLine.posinterp[2,:],mdbm_solutionLine.posinterp[3,:], color="red", leg=false, xlim=(-5,5), ylim=(-5,5), zlim=(-5,5), markersize=2)

@time DTconnect!(mdbm_solutionLine)# connection for plotting (not needed in case of scatter)
plot!(mdbm_solutionLine.posinterp[1,mdbm_solutionLine.DT1'],mdbm_solutionLine.posinterp[2,mdbm_solutionLine.DT1'],mdbm_solutionLine.posinterp[3,mdbm_solutionLine.DT1'], color="green",linewidth=7,leg=false, xlim=(-5,5), ylim=(-5,5), zlim=(-5,5))
