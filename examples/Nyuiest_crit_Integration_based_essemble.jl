using DifferentialEquations
using ForwardDiff
using StaticArrays

using MDBM
using GLMakie

GLMakie.closeall()
GLMakie.activate!(; title="Encircle test")
cmapall = :rainbow

using BenchmarkTools
BenchmarkTools.DEFAULT_PARAMETERS.samples = 100
BenchmarkTools.DEFAULT_PARAMETERS.evals = 1
BenchmarkTools.DEFAULT_PARAMETERS.seconds = 2.50

function D_chareq(p, d, ω)
# 1. The Characteristic Equation
    λ = 1im * ω
    τ = 0.5
    ζ = 0.02

    D = (0.03 * λ^4 + λ^2 + 2 * ζ * λ + 1 + p * exp(-τ * λ) + d * λ * exp(-τ * λ))
    return D
end
const n_power_max = 4
const ω_end = 5000.0 # Max frequency to integrate to (should be high enough for highest power to dominate)

## # Lobe structure - delay turining modell
# function D_chareq(invΩ, w, ω)
#     λ = 1im * ω
#     τ = 2π * invΩ
#     ζ = 0.01
# 
#     H = 1 / (λ^2 + 2 * ζ * λ + 1)
#     D = (1 / H + 1 + w * (1.0 - exp(-τ * λ)))#/ (abs(λ)^2 + 1)
#     return D
# end
# const ω_end = 15000.0 # Max frequency to integrate to (should be high enough for highest power to dominate)
# const n_power_max = 2

# 2. The Core ODE Solver Function
function calculate_unstable_roots_direct(D_func, ω_max, n_power_max, p, d)

    # Define the ODE system: dy/dω = (1/π) * d(arg D)/dω
    #function phase_ode!(dy, y, params, ω)
    global Dminimu=Inf
    function phase_ode_local(y, params, ω)::SVector{1,Float64}
        # Evaluate D(jω) and its derivative dD/dω in a single step using ForwardDiff.Dual
        dual_ω = ForwardDiff.Dual(ω, 1.0)
        res = D_func(p, d, dual_ω)

        # Extract the complex value and its derivative from the Dual result
        D_val = complex(ForwardDiff.value(real(res)), ForwardDiff.value(imag(res)))
        Dminimu=min(Dminimu, abs(D_val))
        dD_dω = complex(ForwardDiff.partials(real(res), 1), ForwardDiff.partials(imag(res), 1))

        # The rate of change of the phase is the imaginary part of (dD/dω) / D
        return SA[imag(dD_dω / D_val)]
    end

    # Initial condition to account for the infinite RHP semicircle
    #y0 =[n_power_max / 2.0]
    y0 = SA[0.0]
    tspan = (0.0, Float64(ω_max))

    # Set up the ODE Problem
    #prob = ODEProblem(phase_ode!, y0, tspan)
    prob = ODEProblem(phase_ode_local, y0, tspan)

    # Solve using a composite solver: 
    # AutoTsit5(Rosenbrock23()) starts with a fast non-stiff solver (Tsit5) 
    # and automatically switches to a stiff solver (Rosenbrock23) when it detects 
    # the sudden sharp gradients near the stability limits.
    # save_everystep=false and saveat=[ω_max] ensures high memory efficiency.
    sol = solve(prob, AutoTsit5(Rosenbrock23()),
        reltol=1e-5,
        abstol=1e-5,
        save_everystep=false,
        saveat=[ω_max])#,dt=1e-4

    # Extract the final integrated value
    Z_raw = -(1.0 / π) * sol.u[end][1] + n_power_max / 2.0

    # Because the number of roots must be an integer, we round the raw integration 
    # result to clean up any tiny numerical tolerances.
    Z_int = round(Int, Z_raw)

    return Z_int, Z_raw
end




# -------------------------------------------------------------------------
# 2. The Core ODE Solver Functions (Skaláris ODE megoldás a gyorsaságért!)
# -------------------------------------------------------------------------


# B) Ensemble verzió vektoros bemenetre (Ezt használjuk a gyors Brute-Force-hoz)
function calculate_unstable_roots_p_vec(D_func, ω_max, n_power_max, params_vec::Vector{Tuple{Float64,Float64}})


    # Define the ODE system: dy/dω = (1/π) * d(arg D)/dω
    #function phase_ode!(dy, y, params, ω)
    function phase_ode_local(y, params, ω)::SVector{1,Float64}
        p, d = params
        # Evaluate D(jω) and its derivative dD/dω in a single step using ForwardDiff.Dual
        dual_ω = ForwardDiff.Dual(ω, 1.0)
        res = D_func(p, d, dual_ω)

        # Extract the complex value and its derivative from the Dual result
        D_val = complex(ForwardDiff.value(real(res)), ForwardDiff.value(imag(res)))
        dD_dω = complex(ForwardDiff.partials(real(res), 1), ForwardDiff.partials(imag(res), 1))

        # The rate of change of the phase is the imaginary part of (dD/dω) / D
        return SA[imag(dD_dω / D_val)]
    end

    # Initial condition to account for the infinite RHP semicircle
    #y0 =[n_power_max / 2.0]
    y0 = SA[0.0]
    tspan = (0.0, Float64(ω_max))

    # Set up the ODE Problem
    #prob = ODEProblem(phase_ode!, y0, tspan)
    base_prob = ODEProblem(phase_ode_local, y0, tspan, params_vec[1])

    # Függvény, ami frissíti az alap problémát a k-adik paraméterpárral
    prob_func(prob, i, repeat) = remake(prob, p=params_vec[i])

    # !! A GYORSÍTÁS KULCSA: output_func !!
    # Az Ensemble szimuláció alapból eltárolná mind a többezer ODESolution objektumot a memóriában.
    # Ez rengeteg GC (Garbage Collection) időt okoz. Ezzel megmondjuk neki, hogy amint kész,
    # csak a végeredményt tartsa meg, a többit azonnal dobja el a memóriából!
    output_func(sol, i) = (sol.u[end][1], false)

    # Ensemble probléma létrehozása
    ensemble_prob = EnsembleProblem(base_prob,
        prob_func=prob_func,
        output_func=output_func)

    # Megoldás EnsembleThreads() használatával
    sim = solve(ensemble_prob, AutoTsit5(Rosenbrock23()), EnsembleThreads(),
        trajectories=length(params_vec),
        batch_size=20, # Kis kötegek a thread-eknek = kevesebb ütemezési overhead
        reltol=1e-5, abstol=1e-5,
        save_everystep=false, saveat=[ω_max])

    # Eredmények kigyűjtése (a sim.u az output_func miatt most már közvetlenül a Float64 értékek vektora)
    Z_raws = -(1.0 / π) * sim.u .+ n_power_max / 2.0
    Z_ints = round.(Int, Z_raws)

    return Z_ints, Z_raws
end

# Test Parameters
p_test = 0.5
d_test = 0.1

println("Pre-compilation")
calculate_unstable_roots_direct(D_chareq, ω_end, n_power_max, p_test, d_test)

# -------------------------------------------------------------------------
# 3. Összehasonlítás: Ciklusonként vs. Ensemble
# -------------------------------------------------------------------------
## ---------------- 3. Brute-Force Stability Chart Execution ----------------



# Pv = LinRange(-2.01, 4.0, 465)
# Dv = LinRange(-2.01, 5.0, 305)
Pv = LinRange(-2.01, 4.0, 90)
Dv = LinRange(-2.01, 5.0, 60)

Z_mat_int = zeros(Int, length(Pv), length(Dv))
Z_mat_raw = zeros(Float64, length(Pv), length(Dv))


println("Calculating stability chart over $(length(Pv)) x $(length(Dv)) grid...direct")
# Using Threads.@threads to parallelize the brute-force sweep
#@benchmark 
@time @inbounds Threads.@threads for i in 1:length(Pv)
    for j in 1:length(Dv)
        Z_int, Z_raw = calculate_unstable_roots_direct(D_chareq, ω_end, n_power_max, Pv[i], Dv[j])
        Z_mat_int[i, j] = max(0, Z_int) # Clamp negative values caused by numerical noise
        Z_mat_raw[i, j] = Z_raw # Clamp negative values caused by numerical noise
    end
end

# --- MÓDSZER 2: Az új Ensemble módszer ---
println("\n2. Futtatás Ensemble szimulációként egyben...")
# Létrehozzuk az összes pontpárt tartalmazó vektort (Tuple vektor)
params_vec = vec([(Pv[i], Dv[j]) for i in 1:length(Pv), j in 1:length(Dv)])

@time Z_ints_vec, Z_raws_vec = calculate_unstable_roots_p_vec(D_chareq, ω_end, n_power_max, params_vec)

# Visszaalakítjuk az 1D eredményvektort 2D mátrixszá
Z_mat_raw_ens = reshape(Z_raws_vec, length(Pv), length(Dv))
Z_mat_int_ens = reshape(Z_ints_vec, length(Pv), length(Dv))
Z_mat_int_ens .= max.(0, Z_mat_int_ens) # Clamp negatív értékek

Z_to_plot = Z_mat_raw_ens

# --- PLOTTOLÁS ---
f = Figure(size=(1700, 600))
ax1 = GLMakie.Axis3(f[1, 1], title="3D View: Unstable Roots", xlabel="p", ylabel="d", zlabel="Z")

surface!(ax1, Pv, Dv, Z_to_plot, colormap=:viridis)

ax2 = GLMakie.Axis(f[1, 2], title="2D Stability Chart (Dark Blue = Stable)", xlabel="p", ylabel="d")
hm = heatmap!(ax2, Pv, Dv, Z_to_plot, colormap=:viridis)
Colorbar(f[1, 2][1, 2], hm, label="Number of Unstable Roots (Z)")

display(f)

# -------------------------------------------------------------------------
# 4. Stable domain curves with MDBM solution
# -------------------------------------------------------------------------


Niter = 4
# Kisebb felbontás a határgörbe kereséséhez
Pv_mdbm = LinRange(-2.01, 4.0, 30)
Dv_mdbm = LinRange(-2.01, 5.0, 20)

println("\n ----------------- Stable domain curves with MDBM solution - direct -----------------")

# MDBM a skaláris fgv hívást fogja használni (A-verzió)
mdbm_wrapper_direct(p, d)::Float64 = calculate_unstable_roots_direct(D_chareq, ω_end, n_power_max, p, d)[2] - 0.5
boundary_mdbm_direct = MDBM_Problem(mdbm_wrapper_direct, [Pv_mdbm, Dv_mdbm])
@time MDBM.solve!(boundary_mdbm_direct, Niter, verbosity=0, interpolationorder=0, checkneighbourNum=1, doThreadprecomp=true, normp=2.0, ncubetolerance=0.7)


xyz_sol = getinterpolatedsolution(boundary_mdbm_direct)

DT1 = MDBM.connect(boundary_mdbm_direct)
edge2plot_xyz = [reduce(hcat, [i_sol[getindex.(DT1, 1)], i_sol[getindex.(DT1, 2)], fill(NaN, length(DT1))])'[:] for i_sol in xyz_sol]
lines!(ax2, edge2plot_xyz..., linewidth=2, color=:red, label="midpoints solution connected")
display(f)




# # Benchmark MDBM
# b = @benchmarkable MDBM.solve!(y, Niter, verbosity=0, interpolationorder=0, checkneighbourNum=1, doThreadprecomp=true, normp=2.0, ncubetolerance=0.7) setup = (y = MDBM_Problem(mdbm_wrapper_direct, [Pv_mdbm, Dv_mdbm]))
# run(b)