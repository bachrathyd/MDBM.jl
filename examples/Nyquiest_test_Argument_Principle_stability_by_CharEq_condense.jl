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
GLMakie.activate!(; title="Argument Principle Stability Analysis with MDBM")
cmapall = :rainbow

include("argument_principle_with_MDBM.jl")

#repalce with your characteristic equation
function D_chareq(p, d, ω)::ComplexF64
    λ = 1im * ω
    τ = 0.5
    ζ = 0.02

    D = (0.03 * λ^4 + λ^2 + 2 * ζ * λ + 1 + p * exp(-τ * λ) + d * λ * exp(-τ * λ)) #/ (abs(λ)^2+ 1)
    return D
end
const n_power_max = 4# it defines the hihest power of the polynomial used to fit the phase curve for the encirclement number calculation.

#coars grid for MDBM evaluation
Pv = LinRange(-2.01, 4.0, 30)
Dv = LinRange(-2.01, 5.0, 20)
omv = LinRange(-0.012, 10.0, 15)
#default coars grid for the argument principle calculation (it is enriched by the solution points of the MDBM)
ω_coars = [omv..., (LinRange(10.01, 10, 10) .^ 2)...]

@time boundary_mdbm, xy_p_uniq, Ncirc, mesh_points, mesh_faces, mesh_colors, edge2plot_xyz = argument_principle_solver_with_MDBM(D_chareq, [Pv, Dv, omv], n_power_max, ω_coars)








# --------------- Plotting ----------------

f = Figure(size=(800, 600))
ax = GLMakie.Axis(f[1, 1], title="Number of unstable roots (N)", xlabel="P", ylabel="D")

# Construct the explicit mesh
m_obj = GeometryBasics.Mesh(mesh_points, mesh_faces)

# Plot with shading disabled to avoid the specular crash
mh = mesh!(ax, m_obj,
    color=round.(mesh_colors),#round.
    colormap=cmapall, colorrange=(0, maximum(Ncirc)),
    shading=NoShading) # Prevents the ComputePipeline error

Colorbar(f[1, 1][1, 2], mh)#, vertical=false)
# Overlay boundaries
lines!(ax, edge2plot_xyz..., color=:white, linewidth=1)
display(f)
