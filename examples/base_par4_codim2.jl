5 + 5
using Revise
using MDBM

using GLMakie



#-----------------------------

function foo_par4_codim2(x, y,z,r)
    x^2.0 + y^2.0 +z^2.0 - r^2.0  ,  x - sin(z*2)
end
# constraint - calculate only the points where the constraint is satisfied (e.g.: on the positiev side)
function c(x, y, z,r)
    y+x - 0.5
end

mymdbm = MDBM_Problem(foo_par4_codim2, [-3.0:2.0:3.0, -3.0:2.0:3.0, -3.0:2.0:3.0, 1.0:2.0], constraint=c)
@time solve!(mymdbm, 4, interpolationorder=1)


f = Figure()

xyzr_sol = getinterpolatedsolution(mymdbm)
scatter(f[1, 1],xyzr_sol[1:3]..., markersize = 6, color = xyzr_sol[4])

# # show the points where the function is evaluated
# xyzr_val = getevaluatedpoints(mymdbm)
# fval=getevaluatedfunctionvalues(mymdbm)
# colors = map(f -> RGBf(
#     sign(f[1])/2 + 0.5,
#     sign(f[2])/2 + 0.5,
#     0.5
# ), fval);
# scatter!(xyzr_val[1:3]...,color=xyzr_val[4], markersize = 3)
# 
# 
# # connecting and plotting the "mindpoints" of the n-cubes
# DT1 = connect(mymdbm)
# edge2plot_xyz = [reduce(hcat, [i_sol[getindex.(DT1, 1)], i_sol[getindex.(DT1, 2)], fill(NaN, length(DT1))])'[:] for i_sol in xyz_sol]
# lines!(edge2plot_xyz..., linewidth=5)
# 
# #TODO: uncomment for the gradinets
#  #plotting the gradintes
# gxyz=getinterpolatedgradient(mymdbm.ncubes,mymdbm)
# arrows!(xyz_sol..., gxyz[1]..., arrowsize = 0.1,lengthscale=0.3,arrowcolor=:blue)
# arrows!(xyz_sol..., gxyz[2]..., arrowsize = 0.1,lengthscale=0.3,arrowcolor=:green)
# 
# 
# 
# 
# 
# #--------------------------- Sub-cube interpolation----------------
# #calcuatin the sub-cubes interpolations stored in the mymdbm.ncubes[i].posinterp
# interpsubcubesolution!(mymdbm)
# #extracting the resutls to from the 
# path2points = extract_paths(mymdbm)
# 
# #extracting the unique points and plotting
# puniq = unique(collect(Iterators.flatten(Iterators.flatten(path2points))))
# scatter(f[1, 2],getindex.(puniq, 1), getindex.(puniq, 2), getindex.(puniq, 3), markersize = 6, color = :green)
# 
# 
# 
# #exctracing the simplexes for each ncube
# flatened_path2points = collect(Iterators.flatten(path2points))
# #eliminating the points with less than 2 points (caused by fininte precision)
# true_truflatened_path2points = flatened_path2points[length.(flatened_path2points).==2]
# 
# lines2plot = [(Point3f(ploc[1]) , Point3f(ploc[2])) for ploc in true_truflatened_path2points]
# linesegments!(lines2plot, linewidth=5)
# 
# 
# 
# f = Figure(resolution = (1000, 600))
# lines(f[1, 1],edge2plot_xyz..., linewidth=5)#mindpoints of the n-cubes
# linesegments(f[1, 2],lines2plot, linewidth=5)#midpoints of the "faces" of the n-cubes
# 
# 
# 