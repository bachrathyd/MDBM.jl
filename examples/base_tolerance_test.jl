
using Revise
using MDBM
using GLMakie

GLMakie.closeall()
GLMakie.activate!(;title = "tolerance-test-: 2 parameters, codimension 1")


f = Figure(size = (800, 600))

## -----------------------------

function foo_par2_codim1(x, y)
    ((x^2.0 + y) - 1.0^2.0)*x #TODO: test with this function, too
end

# MDBM accept points outside the n-cube. In this case there the continuitiy of the manifold is more likely to maintinad.
# However it leads to "parallel solution", that is a neightbouring coube "exptrapolates" properly the solution, thus it is accepted.
# The extra accepted distance is set by ncubetolerance. So the points in a n-cube [-1,1]x[-1,1]x...x[-1,1] at distance (1+ ncubetolerance) is accepted.
#  The distance is coubted be the norm "normp".
# Note, that larger "normp" and smaller "ncubetolerance" leads to clearer results, but it can create discontinutiy in the solution or even lost segments.
# Note, the oppisite leads to a wilder evaluation range around the bracketing ncubes.


normp_vals=[5,100,5000,Inf]
ncubetolerance_vals = [0.005,0.1,2,10]
for (fig_i,(normp, ncubetolerance)) in enumerate(zip(normp_vals, ncubetolerance_vals))
    println("normp: ", normp, " ncubetolerance: ", ncubetolerance)
    
mymdbm = MDBM_Problem(foo_par2_codim1, [-3.1:3.0, -3.1:3.0])

    
ax1=GLMakie.Axis(f[fig_i, 1])
ax2=GLMakie.Axis(f[fig_i, 2])
@time solve!(mymdbm, 4, interpolationorder=1, normp = normp, ncubetolerance = ncubetolerance)#number of refinements - increase it slightly to see smoother results 

empty!(ax1)
# n-cube interpolation
#interpolate!(mymdbm, interpolationorder=1, normp = normp, ncubetolerance = ncubetolerance) 
xy_sol = getinterpolatedsolution(mymdbm)
scatter!(ax1,xy_sol..., markersize = 15, color = :red)

# show the points where the function is evaluated
xy_val = getevaluatedpoints(mymdbm)
fval=getevaluatedfunctionvalues(mymdbm)
scatter!(ax1,xy_val...,color=sign.(fval))

# connecting and plotting the "mindpoints" of the n-cubes
DT1 = connect(mymdbm)
edge2plot_xyz = [reduce(hcat, [i_sol[getindex.(DT1, 1)], i_sol[getindex.(DT1, 2)], fill(NaN, length(DT1))])'[:] for i_sol in xy_sol]
lines!(ax1,edge2plot_xyz..., linewidth=5)




# --------------------------- Sub-cube interpolation----------------

empty!(ax2)
#calcuatin the sub-cubes interpolations stored in the mymdbm.ncubes[i].posinterp
interpsubcubesolution!(mymdbm, normp = normp, ncubetolerance = ncubetolerance)
#extracting the resutls to from the 
path2points = extract_paths(mymdbm)

#extracting the unique points and plotting
puniq = unique(collect(Iterators.flatten(Iterators.flatten(path2points))))
scatter!(ax2,getindex.(puniq, 1), getindex.(puniq, 2))



#exctracing the simplexes for each ncube
flatened_path2points = collect(Iterators.flatten(path2points))
#eliminating the points with less than 2 points (caused by fininte precision)
true_truflatened_path2points = flatened_path2points[length.(flatened_path2points) .== 2]
#plotting the lines between the points
lines2plot = [(Point2f(ploc[1]) , Point2f(ploc[2])) for ploc in true_truflatened_path2points]
linesegments!(ax2,lines2plot)

end
display(f)
#save("numeric_tolerance.png", f, px_per_unit = 2)  