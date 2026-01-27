using MDBM
using GLMakie
using LinearAlgebra
GLMakie.closeall()
GLMakie.activate!(;title = "2 parameters, codimension 1")
#-----------------------------

function foo_par2_codim1(x, y)
    norm([x,y],10.39999)-1.0#((x^2.0 + y) - 1.0^2.0)*x #TODO: test with this function, too
    #y-sin(x*5.0)/(x*5.0)
    #y-sin(1/x)
    #((x^2.0 + y) - 1.0^2.0)*x #TODO: test with this function, too
end

# # refines all bracketing n-cubes
# mymdbm = MDBM_Problem(foo_par2_codim1, [-3.1:3.0, -3.1:3.0])
# @time solve!(mymdbm, 4,checkneighbourNum=0)#number of refinements - increase it slightly to see smoother results 

# refines only n-cubes where the error is greater than 50% betweenthe worst and best error
mymdbm = MDBM_Problem(foo_par2_codim1, [-3.1:3.0, -3.1:3.0])
@time solve!(mymdbm, 18,refinementratio=0.5)#7)#number of refinements - increase it slightly to see smoother results 
#,checkneighbourNum=2
# # refines only n-cubes where the error is greater than an absolute tolerance
# # stops if the number of refinement is reached or all n-cubes are below the tolerance
# mymdbm = MDBM_Problem(foo_par2_codim1, [-3.1:3.0, -3.1:3.0])
# @time solve!(mymdbm, 20,abstol=1e-3)#7)#number of refinements - increase it slightly to see smoother results 


# Plotting the results
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

f


#--------------------------- Sub-cube interpolation----------------
# valid also for differnet n-cube sizes
ax2=GLMakie.Axis(f[1, 2]; xminorticks = mymdbm.axes[1].ticks, yminorticks  = mymdbm.axes[2].ticks, kwargs...)

#calcuatin the sub-cubes interpolations stored in the mymdbm.ncubes[i].posinterp
@time interpsubcubesolution!(mymdbm, normp=50.0, ncubetolerance=0.1)
#extracting the resutls to from the 
path2points = extract_paths(mymdbm)

# #extracting the unique points and plotting
 puniq = unique(collect(Iterators.flatten(Iterators.flatten(path2points))))
 scatter!(getindex.(puniq, 1), getindex.(puniq, 2),label = "subface - solution")



#exctracing the simplexes for each ncube
flatened_path2points = collect(Iterators.flatten(path2points))
#eliminating the points with less than 2 points (caused by fininte precision)
true_truflatened_path2points = flatened_path2points[length.(flatened_path2points) .== 2]
#plotting the lines between the points
lines2plot = [(Point2f(ploc[1]) , Point2f(ploc[2])) for ploc in true_truflatened_path2points]
linesegments!(lines2plot,label = "subface - connection")

display( f)