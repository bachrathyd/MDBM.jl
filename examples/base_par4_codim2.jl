using MDBM
using GLMakie
GLMakie.closeall()
GLMakie.activate!(; title="4 parameters, codimension 2")


#-----------------------------

function foo_par4_codim2(x, y,z,r)
    x^2.0 + y^2.0 +z^2.0 - r^2.0  ,  x - sin(z*2)
end
# constraint - calculate only the points where the constraint is satisfied (e.g.: on the positiev side)
function c(x, y, z,r)
    y+x - 0.5
end

mymdbm = MDBM_Problem(foo_par4_codim2, [-3.0:2.0:3.0, -3.0:2.0:3.0, -3.0:2.0:3.0, 1.0:2.0], constraint=c)
solve!(mymdbm, 4, interpolationorder=1)


f = Figure()

xyzr_sol = getinterpolatedsolution(mymdbm)
scatter(f[1, 1],xyzr_sol[1:3]..., markersize = 6, color = xyzr_sol[4])

display(GLMakie.Screen(), f)
 #--------------------------- Sub-cube interpolation----------------
 ax2 = GLMakie.Axis3(f[1, 2])
 #calcuatin the sub-cubes interpolations stored in the mymdbm.ncubes[i].posinterp
 
 interpsubcubesolution!(mymdbm)
 #extracting the resutls to from the 
 path2points = extract_paths(mymdbm);
 
 #extracting the unique points and plotting
 puniq = unique(collect(Iterators.flatten(Iterators.flatten(path2points))));
  # plotting the surfce along the first 3 dimension
 # scatter!(ax2, getindex.(puniq, 1), getindex.(puniq, 2), getindex.(puniq, 3), markersize=10, color=:green, label="subface - solution")
 
 
 
 #exctracing the simplexes for each ncube
 flatened_path2points = collect(Iterators.flatten(path2points))
 #eliminating the points with less than 2 points (caused by fininte precision)
 true_truflatened_path2points = flatened_path2points[length.(flatened_path2points).==3]
 
 
 #plotting the lines between the points
 n_faces = reshape(1:(3*length(true_truflatened_path2points)), (3, length(true_truflatened_path2points)))'
 vertices_mat = hcat(Iterators.flatten(true_truflatened_path2points)...)
 # plotting the surfce along the first 3 dimension
 mesh!(vertices_mat[1:3,:], n_faces, alpha=0.5, label="subface - local simplex")

 