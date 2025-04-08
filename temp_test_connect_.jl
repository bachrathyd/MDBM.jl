5+5
using Revise
using MDBM

using Plots
theme(:dark)#:vibrant:dracula:rose_pine
plotly()
#const PLOTS_DEFAULTS = Dict(:theme => :wong, :fontfamily => "Computer Modern", :label => nothing, :dpi => 600 )const PLOTS_DEFAULTS = Dict(:theme => :wong, :fontfamily => "Computer Modern", :label => nothing, :dpi => 600 )
default(size=(1500, 1200), titlefont=(15, "times"), legendfontsize=13, guidefont=(12, :white), tickfont=(12, :orange), guide="x", framestyle=:zerolines, yminorgrid=true, fontfamily="Computer Modern", label=nothing, dpi=600)



#Solution of an uncertain implicit equation

function foo3(x,y,z)
    x^2.0+y^2.0+z-2.0^2.0#,x-sin(z)
end
function foo2(x,y)
    x^2.0+y^2.0-2.0^2.0
end

mymdbm=MDBM_Problem(foo3,[-3.0:3.0,-3.0:3.0,-3.0:3.0])
#mymdbm=MDBM_Problem(foo2,[-3.0:3.0,-3.0:3.0])
solve!(mymdbm,0)



interpolate!(mymdbm,interpolationorder=0)
interpolate!(mymdbm,interpolationorder=1)

x_sol,y_sol,z_sol=
aaavyxcv=getinterpolatedsolution(mymdbm)
scatter(x_sol,y_sol,z_sol,ms=1)

#x_sol,y_sol=getinterpolatedsolution(mymdbm)
#scatter(x_sol,y_sol,ms=2)

plot!(xticks =mymdbm.axes[1].ticks,yticks =mymdbm.axes[2].ticks)


DT1=connect(mymdbm)

for edges in DT1
plot!(
    x_sol[[edges...]],
    y_sol[[edges...]],
    z_sol[[edges...]],
    lw=2)
end
plot!()
