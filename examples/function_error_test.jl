using Revise
using MDBM
using GLMakie
using Makie: Relative
using NonlinearSolve

GLMakie.closeall()
GLMakie.activate!(;title = "Error testing - 2 parameters, codimension 1")
#Checking the error 

function foo_par2_codim1(x, y)
    x^4.0 + y^3.0 - 2.0^2.0+x*y
    #((x^2.0 + y) - 1.0^2.0)*x #TODO: test with this function, too
    abs(sin(x*5))-y+x*y
end

 # refines only n-cubes where the error is greater than 50% betweenthe worst and best error
 mymdbm = MDBM_Problem(foo_par2_codim1, [-3.1:3.0, -3.1:3.0])
 @time MDBM.solve!(mymdbm, 20,refinementratio=0.5,global_max_diff_level=10)#7)#number of refinements - increase it slightly to see smoother results 

 # refines only n-cubes where the error is greater than an absolute tolerance
 # stops if the number of refinement is reached or all n-cubes are below the tolerance
 mymdbm = MDBM_Problem(foo_par2_codim1, [-3.1:3.0, -3.1:3.0])
 @time MDBM.solve!(mymdbm, 20,abstol=1e-2,global_max_diff_level=10)#7)#number of refinements - increase it slightly to see smoother results 


 
 @time MDBM.solve!(mymdbm, 2,abstol=1e-2)#7)#number of refinements - increase it slightly to see smoother results 
##

# Plotting the results
f = Figure()
rowsize!(f.layout, 1, Relative(2/3))
#show the final resolution of the grid based on the minorticks
kwargs = (; xminorticksvisible = true, xminorgridvisible = true, yminorticksvisible = true, yminorgridvisible = true)
ax1=GLMakie.Axis(f[1, 1]; xminorticks = mymdbm.axes[1].ticks, yminorticks  = mymdbm.axes[2].ticks, limits = ((-3.1, 3.0), (-3.1, 3.0)), kwargs...)
ax_hist_pos= GLMakie.Axis(f[1, 2], title="Function error", xlabel="‖error‖", ylabel="count", yscale=log10, xscale=log10)
  ax_hist_fun = GLMakie.Axis(f[2, 2], title="Position error", xlabel="‖error‖", ylabel="count", yscale=log10, xscale=log10)


# show the points where the function is evaluated
xy_val = getevaluatedpoints(mymdbm)
fval=getevaluatedfunctionvalues(mymdbm)
s_eval = scatter!(ax1,xy_val...,color=sign.(fval),label = "evaluated", markersize=3,)
Colorbar(f[2, 1], s_eval, label="sign(f)")

# n-cube interpolation
xy_sol = getinterpolatedsolution(mymdbm)
ncsize=[sum(nc.size) for nc in mymdbm.ncubes]
scatter!(ax1,xy_sol...,  color = ncsize,label = "solution")#, markersize=10,marker ='o')

# Error of the function values at the interpolated points
f_mdbroot=foo_par2_codim1.(xy_sol...) # should be close to zero
let data = abs.(f_mdbroot)
    non_zero_data = data[data .> 0]
    if !isempty(non_zero_data)
        min_val = minimum(non_zero_data)
        max_val = maximum(non_zero_data)
        bins = 10 .^ range(log10(min_val), stop=log10(max_val), length=21)
        hist!(ax_hist_pos, data, bins=bins)
    else
        hist!(ax_hist_pos, data)
    end
end
# Output: 1.4142135623730951

# Checking the pointwise error
pos_error = similar(xy_sol[1])
for (i,(x,y)) in enumerate(xy_sol)
    #zip(xy_sol...)
    #println("x: $x , y: $y , at i: $i")
    root = closest_root_2d(foo_par2_codim1, [x, y]) 
    pos_error[i]=norm(root-[x,y])
    #println("Pointwise error: $(pos_error[i])")
end
let data = abs.(pos_error)
    non_zero_data = data[data .> 0]
    if !isempty(non_zero_data)
        min_val = minimum(non_zero_data)
        max_val = maximum(non_zero_data)
        bins = 10 .^ range(log10(min_val), stop=log10(max_val), length=21)
        hist!(ax_hist_fun, data, bins=bins)
    else
        hist!(ax_hist_fun, data)
    end
end

using LinearAlgebra
using ForwardDiff

"""
Finds the root of f(x, y) = 0 closest to the initial guess.
Works for underdetermined systems using the Moore-Penrose pseudo-inverse.
"""
function closest_root_2d(f, guess::Vector{FT}; tol=1e-10, max_iter=20) where FT
    x = copy(guess)
    
    # We wrap the function to take a vector for ForwardDiff compatibility
    f_vec(v) = f(v[1], v[2])
    
    for i in 1:max_iter
        val = f_vec(x)
        
        # Check convergence
        if abs(val) < tol
            return x
        end
        
        # Calculate Gradient (1x2 row vector)
        G = ForwardDiff.gradient(f_vec, x)' 
        
        # Check for singularity (Zero Gradient)
        # Using the logic from our previous discussion
        gnorm = norm(G)
        if gnorm < 1000 * eps(FT)
            # println("Singular gradient at $x. Using fallback.")
            return x .+ FT(1000.0) 
        end
        
        # Newton Step: x_{n+1} = x_n - G⁺ * f(x_n)
        # For a row vector, G \ -val is the same as (G' / norm(G)^2) * -val
        x += G \ -val
    end
    
    return x
end

# --- Usage Example ---

result = closest_root_2d(foo_par2_codim1, [2.0, -0.5]  )
println("Root found at: ", result)
println("Function value at root: ", foo_par2_codim1(result...))