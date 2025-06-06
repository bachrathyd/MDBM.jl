using MDBM
using GLMakie

GLMakie.closeall()
GLMakie.activate!(;title = "2 parameters, codimension 1")
#-----------------------------

function foo_par2_codim1(x, y)
    x^4.0 + y^3.0 - 2.0^2.0
    #((x^2.0 + y) - 1.0^2.0)*x #TODO: test with this function, too
end
mymdbm = MDBM_Problem(foo_par2_codim1, [-3.1:3.0, -3.1:3.0])
@time solve!(mymdbm, 5)#number of refinements - increase it slightly to see smoother results 

f = Figure()
#show the final resolution of the grid based on the minorticks
kwargs = (; xminorticksvisible = true, xminorgridvisible = true, yminorticksvisible = true, yminorgridvisible = true)
ax1=GLMakie.Axis(f[1, 1]; xminorticks = mymdbm.axes[1].ticks, yminorticks  = mymdbm.axes[2].ticks, kwargs...)

# n-cube interpolation
xy_sol = getinterpolatedsolution(mymdbm)
scatter!(xy_sol..., markersize = 15, color = :red,marker ='x',strokewidth=3,label = "solution")

# show the points where the function is evaluated
xy_val = getevaluatedpoints(mymdbm)
fval=getevaluatedfunctionvalues(mymdbm)
scatter!(xy_val...,color=sign.(fval),label = "evaluated")

# connecting and plotting the "mindpoints" of the n-cubes
DT1 = connect(mymdbm)
edge2plot_xyz = [reduce(hcat, [i_sol[getindex.(DT1, 1)], i_sol[getindex.(DT1, 2)], fill(NaN, length(DT1))])'[:] for i_sol in xy_sol]
lines!(edge2plot_xyz..., linewidth=5,label = "midpoints solution connected")


#plotting the gradintes
gxyz=getinterpolatedgradient(mymdbm.ncubes,mymdbm)
arrows!(xy_sol..., gxyz[1]..., arrowsize = 0.01, lengthscale = 0.1,label = "gradient")#    arrowcolor = strength, linecolor = strength)



#--------------------------- Sub-cube interpolation----------------
ax2=GLMakie.Axis(f[1, 2]; xminorticks = mymdbm.axes[1].ticks, yminorticks  = mymdbm.axes[2].ticks, kwargs...)

#calcuatin the sub-cubes interpolations stored in the mymdbm.ncubes[i].posinterp
@time interpsubcubesolution!(mymdbm)
#extracting the resutls to from the 
path2points = extract_paths(mymdbm)

#extracting the unique points and plotting
puniq = unique(collect(Iterators.flatten(Iterators.flatten(path2points))))
scatter!(getindex.(puniq, 1), getindex.(puniq, 2),label = "subface - solution")



#exctracing the simplexes for each ncube
flatened_path2points = collect(Iterators.flatten(path2points))
#eliminating the points with less than 2 points (caused by fininte precision)
true_truflatened_path2points = flatened_path2points[length.(flatened_path2points) .== 2]
#plotting the lines between the points
lines2plot = [(Point2f(ploc[1]) , Point2f(ploc[2])) for ploc in true_truflatened_path2points]
linesegments!(lines2plot,label = "subface - connection")


display(GLMakie.Screen(), f)

axislegend(ax1)
axislegend(ax2)