# ==============================================================================
#  STABILITY ANALYSIS OF DELAYED DYNAMICAL SYSTEMS WITH EXPERIMENTAL FRF DATA
# ==============================================================================
#
#  OBJECTIVE:
#    Computes and visualizes stability charts for a closed-loop
#    dynamical system consisting the characteristic equation D(par1,par2,...,ω).
#    λ = iω (Laplace variable on the imaginary axis)
#    Assume the function in the form, e.g.: for a P-D controller:
#    function D_chareq(P, D, ω)::ComplexF64
#     ...
#    return D::ComplexF64
#
#  ALGORITHM (2-Step):
#    1. Boundary Tracing (MDBM):
#       Uses the Multi-Dimensional Bisection Method (MDBM) to find the 
#       hyper-surface in parameter space (p..., ω) where roots cross the 
#       imaginary axis (Marginal Stability):
#       Re(D(iω)) = 0  AND  Im(D(iω)) = 0
#
#    2. Stability Verification (Argument Principle):
#       Determines the number of unstable roots (N) in each region by computing 
#       the winding number (Cauchy Argument Principle) of the characteristic 
#       function D(iω) over a frequency sweep:
#       N = P + (1/2π) * Δarg(D(iω))
#
#   The key, theat the Argument Principle is computed in a coars grid every where, where the MDBM 
#   evluated the function (mostly around the boundary). To calculate the Argument Principle is  only
#    a coars grid along ω is used, enriched by the soltuions of the setp-1. This way the in the ω axis,
#    contains higher resolution point where it is cruetaial to determine the encirclement number.

#  OUTPUT:
#    - 3D visualization of the stability manifold.
#    - 2D Stability Charts mapping stable (N=0) vs unstable regions. More percisely, the number of unstable roots (N) is plotted as a color map.
#    - Interpolated heatmaps using Delaunay Triangulation.

# ==============================================================================
using MDBM
using GLMakie

GLMakie.closeall()
GLMakie.activate!(; title="Encircle test")
cmapall = :rainbow

include("argument_principle_with_MDBM.jl")

function D_chareq(p, d, ω)::ComplexF64
    λ = 1im * ω
    τ = 0.5
    ζ = 0.02

    D = (0.03 * λ^4 + λ^2 + 2 * ζ * λ + 1 + p * exp(-τ * λ) + d * λ * exp(-τ * λ)) #/ (abs(λ)^2+ 1)
    return D
end
const n_power_max = 4# it defines the hihest power of the polynomial used to fit the phase curve for the encirclement number calculation.

# Step - 1

function mdbm_chareq_wrapper(p...)::Tuple{Float64,Float64}
    D = D_chareq(p...,)
    return real(D), imag(D)
end

## ------------------------ Bifurcation curves with MDBM solution-------------------------
println(" ----------------- Bifurcation curves with MDBM solution -----------------")
Pv = LinRange(-2.01, 4.0, 30)
Dv = LinRange(-2.01, 5.0, 20)
omv = LinRange(-0.012, 10.0, 35)

boundary_mdbm = MDBM_Problem(mdbm_chareq_wrapper, [Pv, Dv, omv])

#@time solve!(boundary_mdbm, 3, verbosity=1) #number of refinements - increase it slightly to see smoother results 
Niter = 4
@time MDBM.solve!(boundary_mdbm, Niter, verbosity=0, checkneighbourNum=2, doThreadprecomp=true, normp=2.0, ncubetolerance=0.3)

#@show boundary_mdbm

f = Figure(size=(1700, 600))
ax1 = GLMakie.Axis3(f[1, 1])
# # n-cube interpolation
xyz_sol = getinterpolatedsolution(boundary_mdbm)
# scatter!(ax1, xyz_sol..., markersize=6, color=:red, marker='x', strokewidth=3, label="solution")

# connecting and plotting the "mindpoints" of the n-cubes
DT1 = MDBM.connect(boundary_mdbm)
edge2plot_xyz = [reduce(hcat, [i_sol[getindex.(DT1, 1)], i_sol[getindex.(DT1, 2)], fill(NaN, length(DT1))])'[:] for i_sol in xyz_sol]
lines!(ax1, edge2plot_xyz..., linewidth=2, label="midpoints solution connected")
display(f)
# Step - 2
## ------------MDBM - updated Encircle number ------------------------
println(" ----------------- MDBM - updated Encircle number -----------------")
# The coars grid for omega values which must be included in any Encirclement calculation
# Higher ω values are added with a quadratic spacing to capture the behavior at high frequencies without excessive points.
#It is neccessary to determine a correnct "close-to-integer" encirclement number, which is crucial for stability analysis.
#Note, the same axis is used as the coars grid of the MDBM.
ω_coars = [omv..., (LinRange(10.01, 10, 10) .^ 2)...]
xy_p_uniq, Ncirc = argument_principle_with_MDBM(D_chareq,boundary_mdbm, ω_coars,n_power_max)

ax_2D_2 = GLMakie.Axis(f[1, 2])
sf = scatter!(ax_2D_2, getindex.(xy_p_uniq, 1), getindex.(xy_p_uniq, 2),
    color=Ncirc, colormap=cmapall, colorrange=(0, maximum(Ncirc)), markersize=10, label="evaluated")
Colorbar(f[1, 3][1, 2], sf)#, vertical=false)


lines!(ax_2D_2, edge2plot_xyz[1], edge2plot_xyz[2], linewidth=4, label="midpoints solution connected")
display(f)


## step - 3 - triangulation for beautiful high res. visualization of the stability chart
println(" ----------------- triangulation for beautiful high res. visualization of the stability chart -----------------")

@time mesh_points, mesh_faces, mesh_colors,edge2plot_xyz = trinangulation_of_MDBM_results(boundary_mdbm, xy_p_uniq, Ncirc)


# 4. Plotting
ax_5 = GLMakie.Axis(f[1, 3], title="Stability Regions (Ncirc)")

# Construct the explicit mesh
m_obj = GeometryBasics.Mesh(mesh_points, mesh_faces)

# Plot with shading disabled to avoid the specular crash
mh = mesh!(ax_5, m_obj,
    color=round.(mesh_colors),#round.
    colormap=cmapall, colorrange=(0, maximum(Ncirc)),
    shading=NoShading) # Prevents the ComputePipeline error

Colorbar(f[1, 3][1, 2], mh)#, vertical=false)
# Overlay boundaries
lines!(ax_5, edge2plot_xyz..., color=:white, linewidth=1)
display(f)
