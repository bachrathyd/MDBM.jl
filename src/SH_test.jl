# SH test

include("MDBM__types.jl")


ax1=Axis([-5,0,3,5],"a")
ax2=Axis(-3:3.0,"b")

mdbmaxes=[ax1,ax2,ax3]

function f(a,b)
    a^2+b^2-2
end
function c(a,b)
    a-b
end

mymdbm=MDBM_Problem(f,[ax1,ax2],constraint=c)
# mymdbm=MDBM_Problem(f,[ax1,ax2])

interpolate!(mymdbm,interpolationorder=1)

for k=1:5
refine!(mymdbm)
interpolate!(mymdbm,interpolationorder=1)
end

using Plots
gr()

a,b=getinterpolatedpoint(mymdbm)
scatter(a,b)
