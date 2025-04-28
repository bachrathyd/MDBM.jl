5 + 5
using Revise
using MDBM

using GLMakie

#using Plots
#theme(:dark)#:vibrant:dracula:rose_pine
#plotly()
##const PLOTS_DEFAULTS = Dict(:theme => :wong, :fontfamily => "Computer Modern", :label => nothing, :dpi => 600 )const PLOTS_DEFAULTS = Dict(:theme => :wong, :fontfamily => "Computer Modern", :label => nothing, :dpi => 600 )
#default(size=(1500, 1200), titlefont=(15, "times"), legendfontsize=13, guidefont=(12, :white), tickfont=(12, :orange), guide="x", framestyle=:zerolines, yminorgrid=true, fontfamily="Computer Modern", label=nothing, dpi=600)

f = Figure()

#Solution of an uncertain implicit equation

function foo4(x, y, z, r)
    x^2.0 + y^2.0 + z - r^2.0, z - y#,x-sin(z)#
end
function foo3(x, y, z)
    x^2.0 + y^2.0 + z - 2.0^2.0, x - sin(z)#,z-y#
end
function foo2(x, y)
    x^2.0 + y^2.0 - 2.0^2.0
end

mymdbm = MDBM_Problem(foo3, [-3.0:3.0, -3.0:3.0, -3.0:3.0])
# mymdbm=MDBM_Problem(foo2,[-3.0:3.0,-3.0:3.0])
solve!(mymdbm, 1)




#interpolate!(mymdbm,interpolationorder=0)
interpolate!(mymdbm, interpolationorder=1)

xyz_sol = getinterpolatedsolution(mymdbm)
scatter(xyz_sol...)

#plot!(xticks =mymdbm.axes[1].ticks,yticks =mymdbm.axes[2].ticks)


DT1 = connect(mymdbm)
edge2plot_xyz = [reduce(hcat, [i_sol[getindex.(DT1, 1)], i_sol[getindex.(DT1, 2)], fill(NaN, length(DT1))])'[:] for i_sol in xyz_sol]

lines!(edge2plot_xyz...)
#for edges in DT1
#plot!(
#    x_sol[[edges...]],
#    y_sol[[edges...]],
#    z_sol[[edges...]])
#end











#-----------------------------
f = Figure()

GLMakie.Axis(f[1, 1])


function foo2(x, y)

    if all([x == 3.0, y == 3.0])
        println("x,y=$x, $y")
        return 0001.0001
    else
        #(x^-1.0 + y^-1.0 - 2.0^2.0)^1.0
        -(-x^2.0 + y^1.0 - 0.0*2.0^2.0)^1.0
        #(x^2.0 + y^2.0 - 2.0^2.0)^1.0
    end
end
mymdbm = MDBM_Problem(foo2, [[0.101021010, 3.0], 0.0015205: 3.0])
solve!(mymdbm, 5)

#mymdbm=MDBM_Problem(foo2,[-5:5,-3:3])
#solve!(mymdbm,2)


#interpolate!(mymdbm,interpolationorder=0)
interpolate!(mymdbm, interpolationorder=1)

xyz_sol = getinterpolatedsolution(mymdbm)
 #scatter(xyz_sol...)
scatter!(xyz_sol...)

xyz_val = getevaluatedpoints(mymdbm)
scatter!(xyz_val...)

getevaluatedfunctionvalues(mymdbm)





xyz_sol = getinterpolatedsolution(mymdbm) 
gxyz=getinterpolatedgradient(mymdbm.ncubes,mymdbm)


arrows!(xyz_sol..., gxyz[1]..., arrowsize = 10, lengthscale = 0.3)#    arrowcolor = strength, linecolor = strength)












 fi = LinRange(0, 2pi, 5000)
 lines!(2 .* sin.(fi), 2 .* cos.(fi))
# #lines!(xticks =mymdbm.axes[1].ticks,yticks =mymdbm.axes[2].ticks)


DT1 = connect(mymdbm)
edge2plot_xyz = [reduce(hcat, [i_sol[getindex.(DT1, 1)], i_sol[getindex.(DT1, 2)], fill(NaN, length(DT1))])'[:] for i_sol in xyz_sol]

lines!(edge2plot_xyz...)
#for edges in DT1
