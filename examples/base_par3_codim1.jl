using Revise
using MDBM
using GLMakie
#-----------------------------

function foo_par3_codim1(x, y, z)
    x^2.0 + y^2.0 + z^2.0 - 2.0^2.0
end

# constraint - calculate only the points where the constraint is satisfied (e.g.: on the positiev side)
function c(x, y, z)
    x^2.0 + y^2.0 - 0.5^2.0
end

mymdbm = MDBM_Problem(foo_par3_codim1, [-3.0:1.0, -1.0:3.0, -1.0:3.0], constraint=c)
#mymdbm = MDBM_Problem(foo_par3_codim1, [-3.0:1.0, -1.0:3.0, -1.0:3.0]) # without constraint

@time solve!(mymdbm, 2)#number of refinements - increase it slightly to see smoother results 



f = Figure(resolution=(1000, 600))

# n-cube interpolation
interpolate!(mymdbm, interpolationorder=1)
xyz_sol = getinterpolatedsolution(mymdbm)
scatter(f[1, 1], xyz_sol..., markersize=15, color=:red)

# show the points where the function is evaluated
xyz_val = getevaluatedpoints(mymdbm)
fval = getevaluatedfunctionvalues(mymdbm)
scatter!(xyz_val..., color=sign.(fval))

# connecting and plotting the "mindpoints" of the n-cubes
DT1 = connect(mymdbm)
edge2plot_xyz = [reduce(hcat, [i_sol[getindex.(DT1, 1)], i_sol[getindex.(DT1, 2)], fill(NaN, length(DT1))])'[:] for i_sol in xyz_sol]
lines!(edge2plot_xyz..., linewidth=5)

# #plotting the gradintes
gxyz = getinterpolatedgradient(mymdbm.ncubes, mymdbm)
arrows!(xyz_sol..., gxyz[1]..., arrowsize=0.1, lengthscale=0.1)#    arrowcolor = strength, linecolor = strength)





#--------------------------- Sub-cube interpolation----------------
#calcuatin the sub-cubes interpolations stored in the mymdbm.ncubes[i].posinterp
interpsubcubesolution!(mymdbm)
#extracting the resutls to from the 
path2points = extract_paths(mymdbm)

#extracting the unique points and plotting
puniq = unique(collect(Iterators.flatten(Iterators.flatten(path2points))))
scatter(f[1, 2], getindex.(puniq, 1), getindex.(puniq, 2), getindex.(puniq, 3), markersize=10, color=:green)



#exctracing the simplexes for each ncube
flatened_path2points = collect(Iterators.flatten(path2points))
#eliminating the points with less than 2 points (caused by fininte precision)
true_truflatened_path2points = flatened_path2points[length.(flatened_path2points).==3]


#plotting the lines between the points
n_faces = reshape(1:(3*length(true_truflatened_path2points)), (3, length(true_truflatened_path2points)))'
vertices_mat = hcat(Iterators.flatten(true_truflatened_path2points)...)
mesh!(vertices_mat, n_faces, alpha=0.5)

#Final plots with diffent methods
ff = Figure(resolution=(1000, 600))
lines(ff[1, 1], edge2plot_xyz..., linewidth=5)
mesh(ff[1, 2], vertices_mat, n_faces;color = :steelblue,shininess=3000,shading=:FastShading )


#Fun with lights

fig = Figure(size = (600, 600))
ax11 = LScene(fig[1, 1], scenekw = (lights = [DirectionalLight(RGBf(0, 0, 0), Vec3f(-1, 0, 0))],))
ax12 = LScene(fig[1, 2], scenekw = (lights = [DirectionalLight(RGBf(1, 1, 1), Vec3f(-1, 0, 0))],))
lights = [
    DirectionalLight(RGBf(0, 0, 0.7), Vec3f(-1, -1, 0)),
    DirectionalLight(RGBf(0.7, 0.2, 0), Vec3f(-1, 1, -1)),
    DirectionalLight(RGBf(0.7, 0.7, 0.7), Vec3f(1, -1, -1))
]
ax21 = LScene(fig[2, 1], scenekw = (lights = lights,))
ax22 = LScene(fig[2, 2], scenekw = (lights = [DirectionalLight(RGBf(4, 2, 1), Vec3f(0, 0, -1))],))
for ax in (ax11, ax12, ax21, ax22)
    mesh!(ax, vertices_mat, n_faces;color = :steelblue )
end
fig



