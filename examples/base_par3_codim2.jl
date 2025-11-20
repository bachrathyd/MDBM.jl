using MDBM
using GLMakie

GLMakie.closeall()
GLMakie.activate!(; title="3 parameters, codimension 2")
#-----------------------------

function foo_par3_codim2(x, y, z)
    x^2.0 + y^2.0 + z^2.0 - 2.0^2.0, x - sin(z * 2)
end
foo_par3_codim2(1.0, 1.0, 1.0)
mymdbm = MDBM_Problem(foo_par3_codim2, [-3.0:3.0, -3.0:3.0, -3.0:3.0])
#@time solve!(mymdbm, 3, verbosity=1) #number of refinements - increase it slightly to see smoother results 
@time MDBM.solve!(mymdbm, 3, verbosity=1, normp=1, ncubetolerance=10, checkneighbourNum=1, doThreadprecomp=true)
println(mymdbm)

f = Figure(size=(1000, 600))
ax1 = GLMakie.Axis3(f[1, 1])
# n-cube interpolation
xyz_sol = getinterpolatedsolution(mymdbm)
scatter!(ax1, xyz_sol..., markersize=6, color=:red, marker='x', strokewidth=3, label="solution")

# show the points where the function is evaluated
xyz_val = getevaluatedpoints(mymdbm)
fval = getevaluatedfunctionvalues(mymdbm)
colors = map(f -> RGBf(
        sign(f[1]) / 2 + 0.5,
        sign(f[2]) / 2 + 0.5,
        0.5
    ), fval);
scatter!(xyz_val..., color=colors, markersize=3, label="evaluated")


# connecting and plotting the "mindpoints" of the n-cubes
DT1 = connect(mymdbm)
edge2plot_xyz = [reduce(hcat, [i_sol[getindex.(DT1, 1)], i_sol[getindex.(DT1, 2)], fill(NaN, length(DT1))])'[:] for i_sol in xyz_sol]
lines!(edge2plot_xyz..., linewidth=5, label="midpoints solution connected")

# #TODO: uncomment for the gradinets
# #plotting the gradintes
# gxyz = getinterpolatedgradient(mymdbm.ncubes, mymdbm)
# arrows!(xyz_sol..., gxyz[1]..., arrowsize=0.1, lengthscale=0.3, arrowcolor=:blue)
# arrows!(xyz_sol..., gxyz[2]..., arrowsize=0.1, lengthscale=0.3, arrowcolor=:green)



#--------------------------- Sub-cube interpolation----------------
ax2 = GLMakie.Axis3(f[1, 2])
#calcuatin the sub-cubes interpolations stored in the mymdbm.ncubes[i].posinterp
interpsubcubesolution!(mymdbm)
#extracting the resutls to from the 
path2points = extract_paths(mymdbm)

#extracting the unique points and plotting
puniq = unique(collect(Iterators.flatten(Iterators.flatten(path2points))))
scatter!(ax2, getindex.(puniq, 1), getindex.(puniq, 2), getindex.(puniq, 3), markersize=6, color=:green, label="subface - solution")



#exctracing the simplexes for each ncube
flatened_path2points = collect(Iterators.flatten(path2points))
#eliminating the points with less than 2 points (caused by fininte precision)
true_truflatened_path2points = flatened_path2points[length.(flatened_path2points).==2]

lines2plot = [(Point3f(ploc[1]), Point3f(ploc[2])) for ploc in true_truflatened_path2points]
linesegments!(lines2plot, linewidth=5, label="subface - connection")


display(GLMakie.Screen(), f)
axislegend(ax1)
axislegend(ax2)


fig = Figure(size=(1000, 600))
lines(fig[1, 1], edge2plot_xyz..., linewidth=5, label="midpoints solution connected")#mindpoints of the n-cubes
linesegments(fig[1, 2], lines2plot, linewidth=5, label="subface - connection")#midpoints of the "faces" of the n-cubes

display(GLMakie.Screen(), fig)
fig
