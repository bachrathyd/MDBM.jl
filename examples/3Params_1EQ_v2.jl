include((x->x[1:find(t->t=='/', x)[end]])(pwd())*"src/MDBM.jl")
using MDBM
using Plots
plotly()
# plotlyjs()

# -----3 parameter , 1 implicit equation ---------------
mdbm_solutionSurf=mdbm_problem(f1,ax) # initialization only

# run the steps separately
for k in 1:2
    refine!(mdbm_solutionSurf)
    checkncube!(mdbm_solutionSurf)
end
checkneighbour!(mdbm_solutionSurf)
interpolate!(mdbm_solutionSurf) # interpolation at the end

scatter(mdbm_solutionSurf.posinterp[1,:],mdbm_solutionSurf.posinterp[2,:],mdbm_solutionSurf.posinterp[3,:], leg=false, xlim=(-5,5), ylim=(-5,5), zlim=(-5,5), markersize=1)

