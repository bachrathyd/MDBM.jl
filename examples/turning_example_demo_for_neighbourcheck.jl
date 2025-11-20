using MDBM
using GLMakie

# Close any existing figures and activate a fresh GLMakie session
GLMakie.closeall()
GLMakie.activate!(; title="Turning with Regeneration — Characteristic Equation")

#------------------------------------------------------------------------------
# 1) CHARACTERISTIC EQUATION of the Regenerative Turning  Model
#------------------------------------------------------------------------------
"""
    char_eq(Ω, w, θ; ζ=0.01)

Residual of the dimensionless chatter characteristic equation

    λ^2 + 2ζ·λ + 1 + w·(1exp(-λ·τ)) = 0

where:
  • Ω  = dimensionless spindle speed
  • w  = dimensionless axial depth of cut
  • ωc = chatter frequency (rad/s)
  • λ  = i·ωc (the characteristic exponent)
  • τ  = 2π / Ω (dimensionless delay)

Returns the squared magnitude |λ^2 + 2ζ·λ + 1 - w·exp(-λ·τ)|²,
which MDBM drives to zero to locate stability boundaries.
"""
function char_eq(Ω, w, ωc; ζ=0.01)
    #ωc =ω              # chatter frequency
    λ = im * ωc              # characteristic exponent
    τ = 2π / Ω               # delay from spindle speed
    D = λ^2 + 2ζ * λ + 1 + w * (1 - exp(-λ * τ))
    return real(D), imag(D)
end
#------------------------------------------------------------------------------
# 2) SET UP MDBM PROBLEM
#------------------------------------------------------------------------------
# Define MDBM axes:
Ω_axis = MDBM.Axis(LinRange(0.1, 1.5, 10), "Ω")# dimensionless spindle speed
w_axis = MDBM.Axis(LinRange(-1.2, 1.2, 5), "w") # axial depth of cut
ωc_axis = MDBM.Axis(LinRange(0.0, 3.0, 7), "ωc") # chatter frequency

#prob = MDBM_Problem((Ω, w, ωc) -> char_eq(Ω, w, ωc; ζ=0.02), [Ω_axis, w_axis, ωc_axis])
turning_mdbm_prob = MDBM_Problem(char_eq, [Ω_axis, w_axis, ωc_axis])

solve!(turning_mdbm_prob, 4, verbosity=1, checkneighbourNum=2)# check neighbour at every iteration
#solve!(turning_mdbm_prob, 4, verbosity=1,checkneighbourNum=1)# (default) check neighbour only once at the end, in some case it is slower
#solve!(turning_mdbm_prob, 4, verbosity=1,checkneighbourNum=0)# 

#------------------------------------------------------------------------------
# 3) SOLVE AND PLOT THE STABILITY SURFACE
#------------------------------------------------------------------------------
# Perform refinements (e.g., 6) to capture the boundary accurately

# Extract the interpolated solution (3 × N array of [Ω; w; θ])
Ω_sol, w_sol, ωc_sol = getinterpolatedsolution(turning_mdbm_prob)


# 3D scatter of the stability boundary in (Ω, w, ωc)
f = Figure(size=(900, 600))
ax = Axis3(f[1, 1];
    xlabel="Ω (ω_s/ω_n)",
    ylabel="w (axial depth)",
    zlabel="ωc (chatter freq.)",
    title="Stability Boundary — Turning Regeneration")
scatter!(ax, Ω_sol, w_sol, ωc_sol;
    markersize=5, color=:tomato, label="roots")
f


#---------------------------- testing the checkneighbourNum -------------------
# This section tests the effect of different `checkneighbourNum` values on the solution.
f_test = Figure(size=(900, 600))
label_test = ["no neighbour check at all", "(default) check neighbour only once at the end, in some case it is slower", "check neighbour at every iteration"]
for checkneighbourNum_test in [0, 1, 2]
    turning_mdbm_prob_test = MDBM_Problem(char_eq, [Ω_axis, w_axis, ωc_axis])
    time_test = @elapsed solve!(turning_mdbm_prob_test, 4, checkneighbourNum=checkneighbourNum_test)# check neighbour at every iteration

    # Extract the interpolated solution (3 × N array of [Ω; w; θ])
   local Ω_sol, w_sol, ωc_sol = getinterpolatedsolution(turning_mdbm_prob_test)


    # 3D scatter of the stability boundary in (Ω, w, ωc)
    ax_loc = Axis3(f_test[1, checkneighbourNum_test+1];
        xlabel="Ω (ω_s/ω_n)",
        ylabel="w (axial depth)",
        zlabel="ωc (chatter freq.)",
        title=label_test[checkneighbourNum_test+1] * " \n checkneighbourNum_test: $checkneighbourNum_test, time: $(round(time_test, digits=2))s")
    scatter!(ax_loc, Ω_sol, w_sol, ωc_sol;
        markersize=5, color=:tomato, label="roots")

    xyz_val = getevaluatedpoints(turning_mdbm_prob_test)
    fval = getevaluatedfunctionvalues(turning_mdbm_prob_test)

    scatter!(ax_loc,xyz_val..., markersize=3, label="evaluated")

end
f_test