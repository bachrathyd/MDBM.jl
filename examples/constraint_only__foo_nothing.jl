# This is a demo file showing how to use MDBM to find the boundary of a region defined by constraints.
# It demonstrates two approaches:
# 1. Using MDBM with only constraints (no equality constraints).
# 2. Defining the boundary as an implicit function (minimum of constraints).

using MDBM
using GLMakie
using LinearAlgebra

GLMakie.closeall()
GLMakie.activate!(; title="Encircle test")

# Setup the plotting environment
f = Figure(size=(1000, 900))
ax11 = GLMakie.Axis(f[1, 1])
ax12 = GLMakie.Axis(f[1, 2])
ax21 = GLMakie.Axis(f[2, 1])
ax22 = GLMakie.Axis(f[2, 2])

# Define the initial parameter space (grid)
x_coars=-2.0:2.0
y_coars=-1.0:2.0

# Define the constraint function
# It returns a tuple of values. The region of interest is where all values are > 0.
function coo_par2(x, y)
    4.0-x.^4 - (y.+0.5)^2, x + 10 * y - 0.4x^2
end


# ==============================================================================
# Section 1: MDBM with constraints only
# ==============================================================================
# In this approach, we provide a "dummy" function `foo` that returns `nothing`.
# MDBM will detect the boundaries of the constraints `coo_par2`.

# Dummy function (no equality constraint)
foo(p...) = nothing

# Create the MDBM problem
mdbm_const_only = MDBM_Problem(foo, [x_coars, y_coars], constraint=coo_par2)

# Solve the problem
# abstol=1e-3 ensures refinement stops when the n-cubes are small enough.
solve!(mdbm_const_only, 10, abstol=1e-3, local_max_diff_level=0)


# -------- Plotting Section 1 --------

# Get the interpolated solution points
xyz_sol = getinterpolatedsolution(mdbm_const_only)

# Connect the solution points to form lines (for visualization)
# `connectoverlap` is preferred for variable size n-cubes (adaptive refinement).
DT1 = MDBM.connectoverlap(mdbm_const_only)

# Prepare line segments for plotting
edge2plot_xyz = [reduce(hcat, [i_sol[getindex.(DT1, 1)], i_sol[getindex.(DT1, 2)], fill(NaN, length(DT1))])'[:] for i_sol in xyz_sol]

# Triangulate the solution for mesh plotting (if needed)
DT2 = triangulation(DT1)
DT2_mat = vcat(transpose.(collect.(DT2))...,);

# Plot on Axis 1,1
empty!(ax11)
ax11.title = "Constraint Only: Steps"
xy_val = getevaluatedpoints(mdbm_const_only)
scatter!(ax11, xy_val..., markersize=4, color=:gray, label="evaluated")
lines!(ax11, edge2plot_xyz..., linewidth=4, color=:magenta, label="connections")
scatter!(ax11, xyz_sol..., markersize=5, color=:black, label="solution")
axislegend(ax11)

# Plot mesh on Axis 2,1
empty!(ax21)
ax21.title = "Constraint Only: Final Result"
mesh!(ax21, hcat(xyz_sol...), DT2_mat, label="boundary")#, color=xyz_sol[1]) # mesh is typically for surfaces
lines!(ax21, edge2plot_xyz..., linewidth=3, label="boundary")

# ==============================================================================
# Section 2: Direct solution of the boundary
# ==============================================================================
# Here we define the boundary as an implicit equation: min(c1, c2, ...) = 0.
# This effectively finds the boundary of the intersection of the constraints.

# Define the implicit function
foo_c_min(x...,)=minimum([coo_par2(x...)...])

# Create the MDBM problem
mdbm_foo = MDBM_Problem(foo_c_min, [x_coars, y_coars])

# Solve
solve!(mdbm_foo, 10, abstol=1e-3, local_max_diff_level=1)

# -------- Plotting Section 2 --------

empty!(ax12)
ax12.title = "Implicit Function: Steps"
xyz_sol = getinterpolatedsolution(mdbm_foo)

# Connect points
DT1 = MDBM.connectoverlap(mdbm_foo)

# Prepare lines
edge2plot_xyz = [reduce(hcat, [i_sol[getindex.(DT1, 1)], i_sol[getindex.(DT1, 2)], fill(NaN, length(DT1))])'[:] for i_sol in xyz_sol]

# Plot on Axis 1,2
xy_val = getevaluatedpoints(mdbm_foo)
scatter!(ax12, xy_val..., markersize=4, color=:gray, label="evaluated")
lines!(ax12, edge2plot_xyz..., linewidth=4, color=:magenta, label="connections")
scatter!(ax12, xyz_sol..., markersize=5, color=:black, label="solution")
axislegend(ax12)

# Plot lines on Axis 2,2
empty!(ax22)
ax22.title = "Implicit Function: Final Result"
lines!(ax22, edge2plot_xyz..., linewidth=3, label="boundary")
display(f)
