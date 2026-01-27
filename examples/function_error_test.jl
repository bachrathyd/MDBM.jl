5+5
using Revise
using MDBM
using GLMakie

GLMakie.closeall()
GLMakie.activate!(; title="Error testing - 2 parameters, codimension 1")
#Checking the error 

function foo_par2_codim1(x, y)
    x^4.0 + y^3.0 - 2.0^2.0 + x * y
    #((x^2.0 + y) - 1.0^2.0)*x #TODO: test with this function, too
    abs(sin(x * 5 - 0.5)) - y + x * y
    #*(sign(0.5*x+y))
    
   # pow= 1.1
   # abs(x)^pow + abs(y)^pow - 2.0^pow
end



mymdbm = MDBM_Problem(foo_par2_codim1, [-3.1:3.0, -3.1:3.0])

# Plotting the results
f = Figure(size=(1600, 1100))
#show the final resolution of the grid based on the minorticks
kwargs = ()#; xminorticks=mymdbm.axes[1].ticks, yminorticks=mymdbm.axes[2].ticks, xminorticksvisible=true, xminorgridvisible=true, yminorticksvisible=true, yminorgridvisible=true)
ax1 = GLMakie.Axis(f[1, 1]; limits=((-3.1, 3.0), (-3.1, 3.0)), kwargs...)

ax_hist_fun = GLMakie.Axis(f[1, 3], yscale=log10, title="Function log-Error histogram", xlabel="log10(‖error‖)", ylabel="count")

# 
# @time MDBM.solve!(mymdbm, 20, abstol=1, global_max_diff_level=10)#7)#number of refinements - increase it slightly to see smoother results 
#
# for pow = 0:-0.25:-4
#     tol = 10.0^pow
#     println("Solving with absolute tolerance: $tol")
#     # refines only n-cubes where the error is greater than an absolute tolerance
#     # stops if the number of refinement is reached or all n-cubes are below the tolerance
#     mymdbm = MDBM_Problem(foo_par2_codim1, [-3.1:3.0, -3.1:3.0])
#     @time MDBM.solve!(mymdbm, 50, abstol=tol, refinementratio=0.5, global_max_diff_level=10)
#     @show mymdbm
#     ##
#

 #   mymdbm = MDBM_Problem(foo_par2_codim1, [-3.1:3.0, -3.1:3.0])
 #   @time MDBM.solve!(mymdbm, 20, abstol=1e-3, refinementratio=0.7)
## 

# if true#false#true#
  #   mymdbm = MDBM_Problem(foo_par2_codim1, [-3.1:3.0, -3.1:3.0])
  #   @time MDBM.solve!(mymdbm, 12, refinementratio=0.4, global_max_diff_level=8)
# else
     mymdbm = MDBM_Problem(foo_par2_codim1, [-3.1:3.0, -3.1:3.0])
     @time MDBM.solve!(mymdbm, 15,  abstol=1.0e-3,refinementratio=0.4, global_max_diff_level=8)
# end

for Niter = 1:0#10
    # refines only n-cubes where the error is greater than an absolute tolerance
    # stops if the number of refinement is reached or all n-cubes are below the tolerance
    ## mymdbm = MDBM_Problem(foo_par2_codim1, [-3.1:3.0, -3.1:3.0])
    ## @time MDBM.solve!(mymdbm, Niter, refinementratio=0.3, global_max_diff_level=10)
    @time MDBM.solve!(mymdbm, 1, refinementratio=0.5, global_max_diff_level=5)
    @show mymdbm
    ##

end
# show the points where the function is evaluated
xy_val = getevaluatedpoints(mymdbm)
fval = getevaluatedfunctionvalues(mymdbm)
empty!(ax1)
scatter!(ax1, xy_val..., color=sign.(fval), label="evaluated", markersize=3,)



# n-cube interpolation
xy_sol = getinterpolatedsolution(mymdbm)
ncsize = [sum(nc.size) for nc in mymdbm.ncubes]
@show unique(sort(ncsize))

# Error of the function values at the interpolated points
f_logerror_mdbmroot = log10.(abs.(foo_par2_codim1.(xy_sol...))) # should be close to zero
s_eval = scatter!(ax1, xy_sol..., color=log10.(ncsize), label="evaluated", markersize=10,)
#s_eval = scatter!(ax1, xy_sol..., color=f_logerror_mdbmroot, label="evaluated", markersize=10,)
Colorbar(f[1, 2], s_eval, label="log(cubesize)", vertical=true)


hist!(ax_hist_fun, f_logerror_mdbmroot, bins=20,alpha=0.7)
# hist!(ax_hist_fun, ncsize)
@show mymdbm
display(f)
#end5