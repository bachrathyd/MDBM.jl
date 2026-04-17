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




## 1. The Characteristic Equation
function D_chareq(p, d, ω::T, σ=zero(ω)) where T<:Number
    λ = 1im * ω + σ
    τ = 0.5
    ζ = 0.02
    D = (0.03 * λ^4 + λ^2 + 2 * ζ * λ + 1 + p * exp(-τ * λ) + d * λ * exp(-τ * λ))
    return D
end

const ω_end = 5000.0 # Max frequency to integrate to (should be high enough for highest power to dominate)
D_chareq(0.2, 0.1, ω_end)
D_chareq(0.2, 0.1, ω_end, 0.0)


# using StaticArrays
# 
# function test_4th_order_model(u, h, params, t)
#     p, d = params
#     τ = 0.5
#     ζ = 0.02
# 
#     # Késleltetett állapotok lekérése
#     u_delayed = h(params, t - τ)
# 
#     # Állapotváltozók deriváltjai
#     du1 = u[2]
#     du2 = u[3]
#     du3 = u[4]
#     
#     # A legmagasabb (4.) derivált kifejezése az egyenletből:
#     # 0.03*x^(4) = -x'' - 2ζx' - x - p*x_tau - d*x'_tau
#     du4 = (1.0 / 0.03) * (-u[3] - 2.0*ζ*u[2] - u[1] - p*u_delayed[1] - d*u_delayed[2])
# 
#     return SA[du1, du2, du3, du4]
# end
# 
# p_test = (0.2, 0.1) # (p, d)
# λ_test = 0.0 + 1im * 5.0 # σ_test + 1im * ω_test
# 
# # Emlékezz: az optimalizált Val{N} miatt Val(4)-et adunk át!
# D_numerical = get_D_from_model(test_4th_order_model, λ_test, p_test, Val(4), verbosity=1)*0.03
# D_chareq(0.2, 0.1,5.0)
# 
# 
# function D_chareq(p, d, ω::T, σ::T=zero(ω)) where T<:Number
#     get_D_from_model(test_4th_order_model, 1im * ω + σ, (p,d), Val(4), verbosity=0)*0.03
# end

D_chareq(0.2, 0.1, 5.0)



# # Lobe structure - delay turining modell
# function D_chareq(invΩ, w, ω)
#     λ = 1im * ω
#     τ = 2π * invΩ
#     ζ = 0.01
#     H = 1 / (λ^2 + 2 * ζ * λ + 1)
#     D = (1 / H + 1 + w * (1.0 - exp(-τ * λ)))#/ (abs(λ)^2 + 1)
#     return D
# end
# const ω_end = 5000.0 # Max frequency to integrate to (should be high enough for highest power to dominate)
# 
# function D_chareq(p, d, ω::T, σ::T=zero(ω)) where T<:Number
#     D_chareq_wrapper(p, d_test, 1im * ω + σ)
#     #D_chareq_analytical(p, d_test, 1im * ω + σ)
# end
# const ω_end = 5000.0 # Max frequency to integrate to (should be high enough for highest power to dominate)
# D_chareq(0.2, 0.1, ω_end)
# D_chareq(0.2, 0.1, ω_end, 0.1)


# using LinearAlgebra
# function D_chareq(p, d, ω)
#     λ = 1im * ω
#     τ = 0.5
#     ζ = 0.02
#     M = @SMatrix [
#         (λ^2+ζ*λ+1.0) (-1.0);
#         (-1.0) (0.1*λ^2+ζ*λ+2.0+p*exp(-τ * λ)+d*λ*exp(-τ * λ))
#     ]
#     # A det() simán megeszi a Dual számokkal teli SMatrixot
#     return det(M)
# end
# const ω_end = 1500.0 # Max frequency to integrate to (should be high enough for highest power to dominate)
# const n_power_max = 4

function get_n_power_max(D_func, p, d; ω_large=1e10, σ=0.0)
    # 1. Beállítunk egy kellően nagy frekvenciát
    T_tag = ForwardDiff.Tag{typeof(D_func),typeof(ω_large)}
    dual_ω = ForwardDiff.Dual{T_tag}(ω_large, 1.0)

    # 2. Kiértékeljük a függvényt kikapcsolt (0, 0) késleltetéssel
    res = D_func(p, d, dual_ω, σ)

    # 3. Kinyerjük az értéket és a deriváltat
    D_val = complex(ForwardDiff.value(real(res)), ForwardDiff.value(imag(res)))
    dD_dω = complex(ForwardDiff.partials(real(res), 1), ForwardDiff.partials(imag(res), 1))

    # 4. Kiszámoljuk a matematikai aszimptotát: ω * Re(D' / D)
    n_est = ω_large * real(dD_dω / D_val)

    # 5. Mivel a fokszám garantáltan egész szám, egyszerűen kerekítjük
    return round(Int, n_est), n_est
end

# 2. The Core ODE Solver Function
function calculate_unstable_roots_direct(D_func, ω_max, p, d, σ=0.0)

    # A legkisebb távolság nyomon követése egy helyi referenciával
    min_D = Ref(Inf)
    estimated_sigma = Ref(Inf) # ÚJ: A kritikus gyök σ értékének becslése

    ω_crit = Ref(0.0)
    # Define the ODE system: dy/dω = (1/π) * d(arg D)/dω
    function phase_ode_local(y, params, ω)::SVector{1,Float64}

        # Explicit Tag generálása a "Nothing" hiba elkerülése végett
        T_tag = ForwardDiff.Tag{typeof(D_func),typeof(ω)}
        # Evaluate D(jω) and its derivative dD/dω in a single step using ForwardDiff.Dual
        dual_ω = ForwardDiff.Dual{T_tag}(ω, 1.0)

        # Mivel a default σ=0, az általad írt új függvény tökéletesen lefut így is!
        res = D_func(p, d, dual_ω, σ)

        # Extract the complex value and its derivative from the Dual result
        D_val = complex(ForwardDiff.value(real(res)), ForwardDiff.value(imag(res)))
        dD_dω = complex(ForwardDiff.partials(real(res), 1), ForwardDiff.partials(imag(res), 1))

        darg = imag(dD_dω / D_val)

        # --- ÚJ: Távolság és gyök-becslés regisztrálása ---


        if abs(D_val) < min_D[]
            min_D[] = abs(D_val)
            ω_crit = ω
            #end
            #if abs(Δσ) < abs(estimated_sigma[])

            # 1. dD/dσ kiszámítása a Cauchy-Riemann egyenletekkel (ingyen megvan!)
            dD_dσ = -1im * dD_dω

            # 2. Newton-Raphson lépés a komplex síkon: Δλ ≈ - D(λ) / D'(λ)
            # Ez megmondja a vektort a képzetes tengelytől a tényleges gyökig
            delta_lambda = -D_val / dD_dσ
            # 3. Mivel jelenleg σ=0-n vagyunk, a gyök valós távolsága maga a Δλ valós része:
            Δσ = real(delta_lambda)

            estimated_sigma[] = Δσ
        end

        # The rate of change of the phase is the imaginary part of (dD/dω) / D
        return SA[darg]
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
    n_power_max = get_n_power_max(D_func, p, d, ω_large=ω_max)[1]
    # Extract the final integrated value
    Z_raw = -(1.0 / π) * sol.u[end][1] + n_power_max / 2.0

    # Because the number of roots must be an integer, we round the raw integration 
    # result to clean up any tiny numerical tolerances.
    Z_int = round(Int, Z_raw)

    return Z_int, Z_raw, min_D[], estimated_sigma[], ω_crit## 
end


# Test Parameters
p_test = 0.5
d_test = 0.1

println("Pre-compilation")
z, z_raw, minD, sigmaest, ω_crit = calculate_unstable_roots_direct(D_chareq, ω_end, p_test, d_test)
D_chareq(p_test, d_test, ω_crit, sigmaest)
# -------------------------------------------------------------------------
# 3. Összehasonlítás: Ciklusonként vs. Ensemble
# -------------------------------------------------------------------------
## ---------------- 3. Brute-Force Stability Chart Execution ----------------



Pv = LinRange(-2.01, 4.0, 90)
Dv = LinRange(-2.01, 5.0, 60)
Pv = LinRange(-2.01, 4.0, 465)
Dv = LinRange(-2.01, 5.0, 305)
# Pv = LinRange(0.025, 0.075, 200)
# Dv = LinRange(1.2, 1.25, 150)

p = 0.445
d = 1.2178
σ_line = -1.061
σ_line = 0.0
Z_mat_int = zeros(Int, length(Pv), length(Dv))
Z_mat_raw = zeros(Float64, length(Pv), length(Dv))
R_mat_raw = zeros(Float64, length(Pv), length(Dv))
σ_mat_raw = zeros(Float64, length(Pv), length(Dv))
ω_mat_raw = zeros(Float64, length(Pv), length(Dv))
Dzero_esti_mat_raw = zeros(Float64, length(Pv), length(Dv))

println("Calculating stability chart over $(length(Pv)) x $(length(Dv)) grid...direct")
# Using Threads.@threads to parallelize the brute-force sweep
#@benchmark 
@time @inbounds Threads.@threads for i in 1:length(Pv)
    @inbounds Threads.@threads for j in 1:length(Dv)
        Z_int, Z_raw, R_min, σ_est, ω_crit = calculate_unstable_roots_direct(D_chareq, ω_end, Pv[i], Dv[j], σ_line)
        Z_mat_int[i, j] = max(0, Z_int) # Clamp negative values caused by numerical noise
        Z_mat_raw[i, j] = Z_raw # Clamp negative values caused by numerical noise
        R_mat_raw[i, j] = minimum([maximum([R_min, -1110.9]), 1110.9]) # Clamp negative values caused by numerical noise
        σ_mat_raw[i, j] = minimum([maximum([σ_est, -0.9]), 1.9])
        ω_mat_raw[i, j] = ω_crit
        Dzero_esti_mat_raw[i, j] = minimum([maximum([abs(D_chareq(Pv[i], Dv[j], ω_crit, σ_est)), -1.0]), 1.0])
    end
end

# Szigmoid telítés, hogy elkerüljük a dimenzió problémát a távoli pontokban
C_scale = 0.5 # Skálázó tényező, szabadon hangolható a meredekséghez

#Z_to_plot =tanh.(C_scale .* R_mat_raw) ./ 2.0
C_to_plot = log.(Z_mat_int)#.+tanh.(C_scale .* R_mat_raw) ./ 2.0
Z_to_plot = log.(Z_mat_int)
C_to_plot = (Z_mat_int)

Z_to_plot = (R_mat_raw)
#C_to_plot=((Z_mat_int == 0) .* 2.0 .- 1.0) .* tanh.(C_scale .* R_mat_raw) ./ 2.0

# --- PLOTTOLÁS ---
f = Figure(size=(1700, 600))
ax1 = GLMakie.Axis3(f[1, 1], title="3D View: Unstable Roots", xlabel="p", ylabel="d", zlabel="Z")

surface!(ax1, Pv, Dv, Z_to_plot, color=C_to_plot, colormap=:viridis)
ax12 = GLMakie.Axis3(f[1, 2], title="3D View: Unstable Roots", xlabel="p", ylabel="d", zlabel="Z")
C_to_plot = 0.0 .* Z_mat_int .+ (Z_mat_int .== 0) .* σ_mat_raw .^ 1
C_to_plot = σ_mat_raw .^ 1
surface!(ax12, Pv, Dv, Z_mat_int, color=C_to_plot, colormap=:viridis)

#surface!(ax1, Pv, Dv, Dzero_esti_mat_raw, color=σ_mat_raw, colormap=:viridis)

ax2 = GLMakie.Axis(f[1, 3], title="2D Stability Chart (Dark Blue = Stable)", xlabel="p", ylabel="d")
hm = heatmap!(ax2, Pv, Dv, C_to_plot, colormap=:viridis)
Colorbar(f[1, 3][1, 2], hm, label="Number of Unstable Roots (Z)")

display(f)
##
# -------------------------------------------------------------------------
# 4. Stable domain curves with MDBM solution
# -------------------------------------------------------------------------
println("\n ----------------- Stable domain curves with MDBM solution - direct -----------------")

@time for σ_loc in LinRange(0, -1.4, 10)
    Niter = 4
    # Kisebb felbontás a határgörbe kereséséhez
    Pv_mdbm = LinRange(-2.01, 4.0, 30)
    Dv_mdbm = LinRange(-2.01, 5.0, 20)



    # MDBM a skaláris fgv hívást fogja használni (A-verzió)
    function mdbm_wrapper_direct(p, d)::Float64
        Z_int, Z_raw, R_min, estimated_sigma, ω_crit = calculate_unstable_roots_direct(D_chareq, ω_end, p, d, σ_loc)

        # Ha a rendszer stabil (Z_int == 0), pozitív a távolság, különben negatív
        sign_val = (max(Z_int, 0) == 0) ? 1.0 : -1.0

        # Szigmoid telítés, hogy elkerüljük a dimenzió problémát a távoli pontokban
        C_scale = 0.5 # Skálázó tényező, szabadon hangolható a meredekséghez
        #   return Z_int .+ (Z_int .== 0) .* estimated_sigma .^ 1
        #  return estimated_sigma
        return sign_val * tanh(C_scale * R_min) / 2.00
        return sign_val * abs(estimated_sigma)
    end

    boundary_mdbm_direct = MDBM_Problem(mdbm_wrapper_direct, [Pv_mdbm, Dv_mdbm])
    @time MDBM.solve!(boundary_mdbm_direct, Niter, verbosity=0, interpolationorder=1, checkneighbourNum=1, doThreadprecomp=true, normp=2.0, ncubetolerance=0.7)


    xyz_sol = getinterpolatedsolution(boundary_mdbm_direct)

    DT1 = MDBM.connect(boundary_mdbm_direct)
    edge2plot_xyz = [reduce(hcat, [i_sol[getindex.(DT1, 1)], i_sol[getindex.(DT1, 2)], fill(NaN, length(DT1))])'[:] for i_sol in xyz_sol]
    lines!(ax2, edge2plot_xyz..., linewidth=2, color=:red, label="midpoints solution connected")
    display(f)

end

display(f)
# # Benchmark MDBM
# b = @benchmarkable MDBM.solve!(y, Niter, verbosity=0, interpolationorder=0, checkneighbourNum=1, doThreadprecomp=true, normp=2.0, ncubetolerance=0.7) setup = (y = MDBM_Problem(mdbm_wrapper_direct, [Pv_mdbm, Dv_mdbm]))
# run(b)

##

# -------------------------------------------------------------------------
# 5. σ-stable domain surfaces with MDBM solution
# -------------------------------------------------------------------------

Niter = 2
# Kisebb felbontás a határgörbe kereséséhez
Pv_mdbm = LinRange(-2.01, 4.0, 30)
Dv_mdbm = LinRange(-2.01, 5.0, 20)
sigma_mdbm = LinRange(-1.01, 1.0, 5)

println("\n ----------------- σ-stable domain surfaces with MDBM solution -----------------")


# MDBM a skaláris fgv hívást fogja használni (A-verzió)
function mdbm_wrapper_direct(p, d, σ)::Float64
    Z_int, Z_raw, R_min = calculate_unstable_roots_direct(D_chareq, ω_end, p, d, σ)

    # Ha a rendszer stabil (Z_int == 0), pozitív a távolság, különben negatív
    sign_val = (max(Z_int, 0) == 0) ? 1.0 : -1.0

    # Szigmoid telítés, hogy elkerüljük a dimenzió problémát a távoli pontokban
    C_scale = 0.5 # Skálázó tényező, szabadon hangolható a meredekséghez

    return sign_val * tanh(C_scale * R_min) / 2.00
end

boundary_mdbm_direct = MDBM_Problem(mdbm_wrapper_direct, [Pv_mdbm, Dv_mdbm, sigma_mdbm])
@time MDBM.solve!(boundary_mdbm_direct, Niter, verbosity=0, interpolationorder=1, checkneighbourNum=1, doThreadprecomp=true, normp=2.0, ncubetolerance=0.7)


xyz_sol = getinterpolatedsolution(boundary_mdbm_direct)

#f = Figure(size=(1700, 600))
#ax4 = GLMakie.Axis3(f[1, 4], title="3D View: Unstable Roots", xlabel="p", ylabel="d", zlabel="Z")
#scatter!(ax12,xyz_sol[1], xyz_sol[2], xyz_sol[3], color=:red, markersize=3, marker='x', strokewidth=2, label="solution")
#scatter!(ax12,xyz_sol[1], xyz_sol[2], -xyz_sol[3], color=:blue, markersize=3, marker='x', strokewidth=2, label="solution")

display(f)


#--------------------------- Sub-cube interpolation----------------
#calcuatin the sub-cubes interpolations stored in the mymdbm.ncubes[i].posinterp
interpsubcubesolution!(boundary_mdbm_direct)
#extracting the resutls to from the 
path2points = extract_paths(boundary_mdbm_direct);

#extracting the unique points and plotting
puniq = unique(collect(Iterators.flatten(Iterators.flatten(path2points))));
#scatter!(ax2, getindex.(puniq, 1), getindex.(puniq, 2), getindex.(puniq, 3), markersize=5, color=:green, label="subface - solution")



#exctracing the simplexes for each ncube
flatened_path2points = collect(Iterators.flatten(path2points))
#eliminating the points with less than 2 points (caused by fininte precision)
true_truflatened_path2points = flatened_path2points[length.(flatened_path2points).==3]


#plotting the lines between the points
n_faces = reshape(1:(3*length(true_truflatened_path2points)), (3, length(true_truflatened_path2points)))'
vertices_mat = hcat(Iterators.flatten(true_truflatened_path2points)...)
mesh!(ax12, vertices_mat, n_faces, alpha=0.5, label="subface - local simplex")



# DT1 = MDBM.connect(boundary_mdbm_direct)
# edge2plot_xyz = [reduce(hcat, [i_sol[getindex.(DT1, 1)], i_sol[getindex.(DT1, 2)], fill(NaN, length(DT1))])'[:] for i_sol in xyz_sol]
# lines!(xyz_sol, edge2plot_xyz..., linewidth=2, color=:red, label="midpoints solution connected")
# lines!(ax12, edge2plot_xyz[1], edge2plot_xyz[2],-1.0 .* edge2plot_xyz[3], linewidth=2, color=:green, label="midpoints solution connected")
display(f)

##
