#
# Compute and plot the (Hopf) bifurcation curves/surfaces of the
# Takács–Stépán stretched‐string tyre‐delay model (Takács et al., EJMS 2009).
#  • First: codimension‐2 curves in (V,L) for fixed Σ=1.8, ζ=0  (Figure 4(a)).
#  • Then: codimension‐2 surface in (V,L,ζ) for fixed Σ=1.8        (lifting by ζ).
#
# Julia 1.8+ MDBM.jl GLMakie.jl
#
using MDBM                        # Multi‐Dimensional Bisection Method 
using GLMakie                     # High‐performance plotting             
using LinearAlgebra               # for complex operations (norm, etc.)

GLMakie.closeall()
GLMakie.activate!(;title = "Delay effects in shimmy dynamics of wheels with stretched string-like tyres - D Takács, G Orosz, G Stépán")
# =============================================================================
# 1. Load DIMENSIONAL PARAMETERS (Table 1 of Takács et al. 2009)
# =============================================================================
# Tyre contact half‐length
const a_dim = 0.04    # [m]
# Tyre relaxation length
const σ_dim = 0.072   # [m]
# Lateral tyre stiffness per unit length
const k_dim = 53506   # [N/m^2]
# Lateral tyre damping per unit length
const b_dim = 140     # [Ns/m^2]
# Undamped natural angular frequency (measured)
const ωn_dim = 15.29   # [rad/s]
# Corresponding natural frequency f_n (Hz)
const f_n = 2.43    # [Hz]   # (not used directly below)

# =============================================================================
# 2. Compute DIMENSIONLESS PARAMETERS from Eq. (21) & (26):
#
#   Σ = σ / a ,     ζ = (ω_n b)/(2 k) ,
#   V = (v/(2a))/ω_n ,    L = l/a .
#
#  We will fix Σ = 1.8 and later vary ζ as needed.
# =============================================================================
const Σ_target = σ_dim / a_dim                # = 1.8  (dimensionless relaxation length)
const ζ_target = (ωn_dim * b_dim) / (2 * k_dim)  # ≈ 0.02 (dimensionless damping ratio)

println()
println("Dimensional → Dimensionless:")
println("  a = $(a_dim) m  ⇒   Σ = σ/a = $(round(Σ_target, digits=3))")
println("  b = $(b_dim) Ns/m^2,  k = $(k_dim) N/m^2,  ωn = $(ωn_dim) rad/s")
println("  ⇒ ζ = ω_n b/(2k) = $(round(ζ_target, digits=4))")
println()

# =============================================================================
# 3. Define the DIMENSIONLESS characteristic function D(λ; μ)  (Eq. 31)
#
#   μ = (V, L, Σ, ζ) 
#   λ = i·ω   (purely imaginary for hopf)
#   D(λ; μ) = Σ V^2 λ^3 
#          + 2 V (V + Σ ζ) λ^2 
#          + (Σ + 4ζV) λ 
#          + 2 
#            − (L − 1 − Σ)/(L^2 + 1/3 + Σ(L^2 + 1 + Σ))
#            × {  2 λ^2 [ (L−1)λ + 2 − ( (L+1)λ + 2 ) e^{−λ} ]
#                + 4ζ V L (1 + Σ) (2 + Σ λ)/(L−1−Σ)
#                + (L − 1 − Σ)(2Σζ V λ + Σ + 4ζ V)
#                + (L + 1 + Σ)(2Σζ V λ + Σ − 4ζ V) e^{−λ}  }
#
#  (All of that is straight from Eq. (31) in Takács et al. (2009).) 
# =============================================================================
"""
    characteristic_D(λ, V, L, Σ, ζ)


# Arguments
- `λ::ComplexF64`: the complex frequency (e.g., `im * ω`).
- `V::Float64`: dimensionless towing-speed parameter.
- `L::Float64`: dimensionless caster-length parameter.
- `Σ::Float64`: dimensionless relaxation‐length (`σ / a`).
- `ζ::Float64`: dimensionless damping ratio (`ω_n b / (2 k)`).

# Returns
- `ComplexF64`: the value of D(λ; μ).
"""
function characteristic_D(λ::ComplexF64, V::Float64, L::Float64, Σ::Float64, ζ::Float64)::ComplexF64
    # Denominator: L^2 + 1/3 + Σ (L^2 + 1 + Σ)
    denom = L^2 + 1 / 3 + Σ * (L^2 + 1 + Σ)
    # Numerator factor outside the curly braces: (L - 1 - Σ)
    num_fac = L - 1 - Σ

    # --- Terms inside the { … } from the provided formula ---
    # Term A: 2/λ^2 * [ (L - 1)λ + 2 - ((L + 1)λ + 2) e^{-λ} ]
    termA = (2.0 / λ^2) * ((L - 1) * λ + 2 - ((L + 1) * λ + 2) * exp(-λ))

    # Term B: 4 ζ V L (1 + Σ) (2 + Σ λ) / (L - 1 - Σ)
    termB = 4.0 * ζ * V * L * (1 + Σ) * (2 + Σ * λ) / num_fac

    # Term C: (L - 1 - Σ)( 2 Σ ζ V λ + Σ + 4 ζ V )
    termC = num_fac * (2 * Σ * ζ * V * λ + Σ + 4 * ζ * V)

    # Term D: (L + 1 + Σ)( 2 Σ ζ V λ + Σ - 4 ζ V ) e^{−λ}
    termD = (L + 1 + Σ) * (2 * Σ * ζ * V * λ + Σ - 4 * ζ * V) * exp(-λ)

    # Combine the curly‐brace terms:
    curly = (termA + termB + termC + termD)

    # --- Polynomial part: Σ V^2 λ^3 + 2 V (V + Σ ζ) λ^2 + (Σ + 4 ζ V) λ + 2 ---
    poly = Σ * V^2 * λ^3
    poly += 2.0 * V * (V + Σ * ζ) * λ^2
    poly += (Σ + 4.0 * ζ * V) * λ
    poly += 2.0

    # Final D(λ; μ)
    return poly - num_fac / denom * curly
end



#  → We will always call it as D(iω; V, L, Σ, ζ).  For convenience:
function char_fun_dimless(V::Float64, L::Float64, ω::Float64, Σ::Float64, ζ::Float64)::Tuple{Float64,Float64}
    λ = im * ω
    Dλ = characteristic_D(λ, V, L, Σ, ζ)
    return real(Dλ), imag(Dλ)
end

# Quick sanity check (must be a Complex):
#   e.g. D(i·1; V=2, L=3, Σ=1.8, ζ=0.02) 
#   just to confirm no errors in transcription:
_begin_sanity = char_fun_dimless(2.0, 3.0, 1.0, Σ_target, ζ_target)
println("Sanity check: D(i·1; V=2, L=3, Σ=1.8, ζ=0.02) = (Re,Im) = $(_begin_sanity)")
println()

# =============================================================================
# 4. CODIMENSION‐2 ROOT‐FINDING in (V,L,ω) with Σ=1.8, ζ=0.0
#
#    Re D(iω; V,L,Σ=1.8, ζ=0) = 0
#    Im D(iω; V,L,Σ=1.8, ζ=0) = 0
#
#  → This yields a 1‐dimensional manifold of solutions in 3D.
#  → We will use MDBM_Problem(char_fun3d, axes3d) to approximate it.
# =============================================================================
# Fix Σ and ζ for Figure 4(a):
const Σ_fixed = 1.8           # = 1.8
const ζ_fixed = 0.000                # ζ = 0 for the undamped case (Fig 4(a))

# Define a small wrapper for MDBM that only uses (V,L,ω):
function char_fun3d(x...)
    V, L, ω = x
    return char_fun_dimless(V, L, ω, Σ_fixed, ζ_fixed)
end

# Choose parameter ranges.  Fig 4(a) typically shows V ∈ [0, 4],  L ∈ [0.5, 5.0].
#   → We follow the same window: 
V_min, V_max = 0.001, 0.6
L_min, L_max = -0.2 + pi / 1000, 7.0
ω_min, ω_max = 0.1, 30.0         # ω usually only goes up to 10 or so in figure

# Build coarse MDBM axes in (V, L, ω):
axis_V3d = V_min:0.02:V_max    # coarse spacing of 0.2 in V
axis_L3d = L_min:0.2:L_max    # coarse spacing of 0.2 in L
axis_ω3d = ω_min:0.4:ω_max    # spacing 0.5 in ω

axes3d = [axis_V3d, axis_L3d, axis_ω3d]

println("→ Setting up MDBM for codim‐2 curve in (V,L,ω) with Σ=1.8, ζ=0 …")
mdbm3d = MDBM_Problem(char_fun3d, axes3d, memoization=true)

# Perform a few refinements (e.g. 4 or 5) to resolve curves:
refine_iters3d = 3

solve!(mdbm3d, refine_iters3d, doThreadprecomp=false, verbosity=1)
println("→ Finished MDBM solve for Fig 4(a).")


# =============================================================================
# 5. PLOTTING CODIM‐2 CURVES in (V,L), colored by ω
#
#   We simply project each triangle (in 3D) down to (V,L), and color by the
#   average ω of its three vertices.  This reproduces the line‐curves in Fig 4(a).
# =============================================================================
println("→ Plotting 2D projection of codim‐2 (V,L) curves, colored by ω …")
GLMakie.activate!()

figure1 = Figure(size=(900, 600))
#ax1 = Axis(figure1[1, 1], xlabel="V (dimensionless)", ylabel="L (dimensionless)", title="Figure 4(a) from Takács et al. (2009), Σ=1.8, ζ=0")

xyzr_sol = getinterpolatedsolution(mdbm3d)
scatter(figure1[1, 1], xyzr_sol[1:2]..., markersize=6, color=xyzr_sol[3])

display(figure1)
save("Fig4a_Takacs_TemporalHopfCurves.png", figure1)

println("→ Saved Figure 4(a) projection as 'Fig4a_Takacs_TemporalHopfCurves.png'.")
println()

# =============================================================================
# 6. LIFTING THE PROBLEM ONE DIMENSION HIGHER:  ADD ζ AS A PARAMETER
#
#   Now treat μ = (V, L, ω, ζ) with Σ fixed at 1.8.  Solve:
#     Re D(iω; V,L,Σ=1.8, ζ) = 0,
#     Im D(iω; V,L,Σ=1.8, ζ) = 0
#   → codim‐2 in 4D ⇒ 2‐dimensional surface in (V,L,ζ,ω).  
#   We will extract that 2D surface via MDBM_Problem(char_fun4d, axes4d).
#
#   Then we project each little triangle from 4D down to 3D (V,L,ζ) and color by ω.
#   This gives a **“surface” in (V,L,ζ)** space, exactly what the user requested.
# =============================================================================
# Build the 4D wrapper char_fun4d:
function char_fun4d(x...)
    V, L, ω, ζ = x
    return char_fun_dimless(V, L, ω, Σ_target, ζ)
end

# Choose ζ range: from 0.0 up to, say, 0.05 to see damping effects (like Fig 4(c)):
ζ_min, ζ_max = 0.00, 0.05

# Define coarse axes in (V, L, ω, ζ):

axis_V4d = V_min:0.1:V_max
axis_L4d = L_min:1.5:L_max
axis_ω4d = ω_min:3.0:ω_max
axis_ζ4d = ζ_min:0.01:ζ_max

axes4d = [axis_V4d, axis_L4d, axis_ω4d, axis_ζ4d]

println("→ Setting up MDBM for codim‐2 SURFACE in (V,L,ω,ζ) …")
mdbm4d = MDBM_Problem(char_fun4d, axes4d)

# Refine a few times:
refine_iters4d = 3   # fewer iters to control runtime
solve!(mdbm4d, refine_iters4d)
println("→ Finished MDBM solve for 4D (V,L,ω,ζ) hopf surface.")


println("→ Plotting 2D projection of codim‐2 (V,L) curves, colored by ω …")



# =============================================================================
# 7. PLOTTING the 3D SURFACE in (V, L, ζ), colored by ω
#
#   This shows, for Σ=1.8, the entire Hopf‐bifurcation SURFACE as ζ varies.
# =============================================================================
println("→ Plotting codim‐2 surface in (V,L,ζ), colored by ω …")
figure2 = Figure(size=(1000, 700))
ax = Axis3(figure2[1, 1],
    xlabel="V (dimless)", ylabel="L (dimless)", zlabel="ζ (dimless)",
    title="Takács Shimmy Hopf Surface: Σ=1.8, 0 ≤ ζ ≤ 0.05",
    aspect=:equal,    # use data‐ratio so the box matches the ranges
    viewmode=:fit      # fit the resulting cuboid into the frame
)
#ax1 = Axis(figure1[1, 1], xlabel="V (dimensionless)", ylabel="L (dimensionless)", title="Figure 4(a) from Takács et al. (2009), Σ=1.8, ζ=0")

xyzr_sol = getinterpolatedsolution(mdbm4d)
scatter!(ax, xyzr_sol[[1, 2, 4]]..., markersize=6, color=xyzr_sol[3])



#--------------------------- Sub-cube interpolation----------------
println("Calculation the interpolated values in the sub-cubes (faces of the n-cubes)")
@time interpsubcubesolution!(mdbm4d)
#extracting the resutls to from the 
path2points = extract_paths(mdbm4d)

# #extracting the unique points and plotting
#puniq = unique(collect(Iterators.flatten(Iterators.flatten(path2points))))
#scatter!(ax, getindex.(puniq, 1), getindex.(puniq, 2), getindex.(puniq, 4), markersize=10, color=getindex.(puniq, 3), label="subface - solution")



#exctracing the simplexes for each ncube
flatened_path2points = collect(Iterators.flatten(path2points))
#eliminating the points with less than 2 points (caused by fininte precision)
true_truflatened_path2points = flatened_path2points[length.(flatened_path2points).==3]


#plotting the lines between the points
n_faces = reshape(1:(3*length(true_truflatened_path2points)), (3, length(true_truflatened_path2points)))'
vertices_mat = hcat(Iterators.flatten(true_truflatened_path2points)...)


ax = Axis3(figure2[1, 2],
    xlabel="V (dimless)", ylabel="L (dimless)", zlabel="ζ (dimless)",
    title="Takács Shimmy Hopf Surface: Σ=1.8, 0 ≤ ζ ≤ 0.05",
    aspect=:equal,    # use data‐ratio so the box matches the ranges
    viewmode=:fit      # fit the resulting cuboid into the frame
)#calcuatin the sub-cubes interpolations stored in the mdbm4d.ncubes[i].posinterp
mesh!(vertices_mat[[1, 2, 4], :], n_faces, alpha=1.0, color=vertices_mat[3, :], label="subface - local simplex")



display(figure2)
save("Takacs_HopfSurface_Σ1p8_VLζ.png", figure2)

println("→ Saved 3D Hopf‐surface as 'Takacs_HopfSurface_Σ1p8_VLζ.png'.")
println()
println("**Done**. You now have:")
println("  • 'Fig4a_Takacs_TemporalHopfCurves.png'   ← Fig 4(a) reconstruction")
println("  • 'Takacs_HopfSurface_Σ1p8_VLζ.png'       ← full 3D surface in (V,L,ζ), colored by ω")
println()



# =============================================================================
# End of script.
# =============================================================================
