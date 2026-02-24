

using MDBM
using SemiDiscretizationMethod

using LinearAlgebra
using StaticArrays

using GLMakie

GLMakie.closeall()
GLMakie.activate!(; title="Stability of Delay Mathieu test")

# Setup the plotting environment
f = Figure(size=(1000, 900))
ax11 = GLMakie.Axis(f[1, 1])


function createMathieuProblem(δ, ε, b0, a1; T=2π)
    AMx = ProportionalMX(t -> @SMatrix [0.0 1.0; -δ-ε*cos(2π / T * t) -a1])
    τ1 = t -> 2π # if function is needed, the use τ1 = t->foo(t)
    BMx1 = DelayMX(τ1, t -> @SMatrix [0.0 0.0; b0 0.0])
    cVec = Additive(t -> @SVector [0.0, sin(4π / T * t)])
    LDDEProblem(AMx, [BMx1], cVec)
end;

τmax = 2π # the largest τ of the system
T = 2π #Principle period of the system (sin(t)=sin(t+P)) 
mathieu_lddep = createMathieuProblem(3.0, 0.2, -0.15, 0.1, T=T); # LDDE problem for Delay Mathieu equation
#--------- Stability map ----------

a1 = 0.1;
ε = 1;
τmax = 2π;
T = 1π;
method = SemiDiscretization(2, T / 40);

foo_δb0(δ, b0)::Float64 = -log(spectralRadiusOfMapping(DiscreteMapping_LR(createMathieuProblem(δ, ε, b0, a1, T=T), method, τmax,
        n_steps=Int((T + 100eps(T)) ÷ method.Δt)), nev=1, tol=1e-3)); # No additive term calculated


axis = [-1:0.5:5.0,
   -2:0.5:1.5]

foo_null(p...) = nothing

DelayMathieu_δb0_mdbm=MDBM_Problem(foo_null, axis,constraint=foo_δb0)# stable area
#DelayMathieu_δb0_mdbm=MDBM_Problem(foo_δb0, axis)# stability border only - almost the same speed
chat_abstol=1e-3
println("Start MDBM solver")
@time solve!(DelayMathieu_δb0_mdbm, 10,abstol=chat_abstol, doThreadprecomp=false,local_max_diff_level=0,verbosity=0)#parallel run must be switched off for eigen valuse calculation
println("Finish MDBM solver")
# -------- Plotting Section 1 --------
# Get the interpolated solution points
xyz_sol = getinterpolatedsolution(DelayMathieu_δb0_mdbm)

# Connect the solution points to form lines (for visualization)
# `connectoverlap` is preferred for variable size n-cubes (adaptive refinement).
DT1 = MDBM.connectoverlap(DelayMathieu_δb0_mdbm)

# Prepare line segments for plotting
edge2plot_xyz = [reduce(hcat, [i_sol[getindex.(DT1, 1)], i_sol[getindex.(DT1, 2)], fill(NaN, length(DT1))])'[:] for i_sol in xyz_sol]

# Triangulate the solution for mesh plotting (if needed)
DT2 = triangulation(DT1)
DT2_mat = vcat(transpose.(collect.(DT2))...,);

# Plot on Axis 1,1
empty!(ax11)

mesh!(ax11, hcat(xyz_sol...), DT2_mat, color=:green)#, color=xyz_sol[1])
lines!(ax11, edge2plot_xyz..., linewidth=1, color=:black, label="midpoints solution connected")
#scatter!(ax11, xyz_sol..., markersize=3, color=:black, label="solution")


# xyz_val = getevaluatedpoints(DelayMathieu_mdbm)
# fval = getevaluatedconstraintvalues(DelayMathieu_mdbm)
# scatter!(ax11,xyz_val..., color=(fval), label="evaluated")

display(f)


#--------- Stability map of the Mathieu equation (no delay)----------

a1 = 0.01;
τmax = 2π;
T = 2π;
b0 = 0.0
method = SemiDiscretization(2, T / 40);

foo_δε(δ, ε)::Float64 = -log(spectralRadiusOfMapping(DiscreteMapping_LR(createMathieuProblem(δ, ε, b0, a1, T=T), method, τmax,
        n_steps=Int((T + 100eps(T)) ÷ method.Δt)), nev=1, tol=1e-4)); # No additive term calculated

axis = [MDBM.Axis(-2:1:5.0, :δ),
    MDBM.Axis(-0.01:1:5, :ε)]


DelayMathieu_δε_mdbm=MDBM_Problem(foo_null, axis,constraint=foo_δε)
coordinate_abstol=1e-2
println("Start MDBM solver")
@time solve!(DelayMathieu_δε_mdbm, 10,abstol=coordinate_abstol, doThreadprecomp=false,local_max_diff_level=0,verbosity=0)#parallel run must be switched off for eigen valuse calculation
println("Finish MDBM solver")
# -------- Plotting Section 1 --------

# Get the interpolated solution points
xyz_sol = getinterpolatedsolution(DelayMathieu_δε_mdbm)

# Connect the solution points to form lines (for visualization)
# `connectoverlap` is preferred for variable size n-cubes (adaptive refinement).
DT1 = MDBM.connectoverlap(DelayMathieu_δε_mdbm)

# Prepare line segments for plotting
edge2plot_xyz = [reduce(hcat, [i_sol[getindex.(DT1, 1)], i_sol[getindex.(DT1, 2)], fill(NaN, length(DT1))])'[:] for i_sol in xyz_sol]

# Triangulate the solution for mesh plotting (if needed)
DT2 = triangulation(DT1)
DT2_mat = vcat(transpose.(collect.(DT2))...,);

# Plot on Axis 1,1

ax12 = GLMakie.Axis(f[1, 2])

mesh!(ax12, hcat(xyz_sol...), DT2_mat, color=:green)#, color=xyz_sol[1])
lines!(ax12, edge2plot_xyz..., linewidth=1, color=:black, label="midpoints solution connected")
#scatter!(ax12, xyz_sol..., markersize=3, color=:black, label="solution")


# xyz_val = getevaluatedpoints(DelayMathieu_mdbm)
# fval = getevaluatedconstraintvalues(DelayMathieu_mdbm)
# scatter!(ax12,xyz_val..., color=(fval), label="evaluated")

display(f)


#--------- 3D stability map ----------
using MDBM
plotly()
a1 = 0.01;
τmax = 2π;
T = 2π;
method = SemiDiscretization(2, T / 40);

foo(δ, b0, ε) = log(spectralRadiusOfMapping(DiscreteMapping_LR(createMathieuProblem(δ, ε, b0, a1, T=T), method, τmax,
        n_steps=Int((T + 100eps(T)) ÷ method.Δt)), nev=1, tol=1e-4)); # No additive term calculated

axis = [Axis(-2:0.5:5.0, :δ),
    Axis(-2:0.5:1.5, :b0),
    Axis(-0.01:0.5:5, :ε)]

iteration = 2;
@time stab_border_points = getinterpolatedsolution(solve!(MDBM_Problem(foo, axis), iteration));

scatter(stab_border_points...,
    label="", title="Stability border of the delay Mathieu equation", xlabel=L"\delta", ylabel=L"b_0", zlabel=L"\epsilon",
    guidefontsize=14, tickfont=font(10), markersize=1, markerstrokewidth=0)

