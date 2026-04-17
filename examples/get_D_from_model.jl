using LinearAlgebra
using ForwardDiff
using StaticArrays


# ==============================================================================
# 1. AZ IDŐTARTOMÁNYBELI DDE MODELL (1. rendű vektoros alak)
# ==============================================================================
function turning_model(u, h, p, t)
    invΩ, w = p
    τ = 2π * invΩ
    ζ = 0.01

    u_delayed = h(p, t - τ)

    du1 = u[2]
    du2 = -(2.0 + w)*u[1] - 2.0*ζ*u[2] + w*u_delayed[1]

    return SA[du1, du2]
end

function D_chareq_analytical(invΩ, w, λ::Complex)
    τ = 2π * invΩ
    ζ = 0.01
    H = 1 / (λ^2 + 2 * ζ * λ + 1)
    D = (1 / H + 1 + w * (1.0 - exp(-τ * λ)))
    return D
end

# ==============================================================================
# 2. A FEKETEDOBOZ KARAKTERISZTIKUS EGYENLET KINYERŐ (STATIC & ALLOCATION-FREE)
# ==============================================================================
# ÚJ: N helyett Val{N}-et használunk, ami lehetővé teszi a tökéletes SMatrix optimalizációt a fordítónak
function get_D_from_model(bc_model, λ::Complex, p, ::Val{N}; u_eq = zeros(SVector{N, Float64}), verbosity::Int = 0) where {N}
    
    # T_val alkalmazkodik a külső ForwardDiff rétegekhez (ha λ duális szám, T_val is az lesz)
    T_val = typeof(real(λ))
    
# --- 1. PERTURBÁCIÓK ELŐÁLLÍTÁSA (SZABVÁNYOS TAG HASZNÁLATÁVAL) ---
    # Automatikusan generálunk egy belső címkét a bc_model alapján
    T_inner_tag = ForwardDiff.Tag{typeof(bc_model), T_val}
    
    A_duals = SVector{N, ForwardDiff.Dual{T_inner_tag, T_val, N}}(
        ntuple(i -> ForwardDiff.Dual{T_inner_tag, T_val, N}(
            zero(T_val), 
            ForwardDiff.Partials{N, T_val}(ntuple(j -> i == j ? one(T_val) : zero(T_val), Val(N)))
        ), Val(N))
    )
    
    # Mivel t=0, exp(0) = 1
    u_current = u_eq .+ A_duals 
    
    # Kamu történeti függvény (h_mock)
    h_mock(p_arg, t_arg; idxs=nothing) = begin
        exp_term = exp(λ * t_arg)
        if idxs === nothing
            return u_eq .+ A_duals .* exp_term
        else
            return u_eq[idxs] + A_duals[idxs] * exp_term
        end
    end
    
    # Modell kiértékelése
    du = bc_model(u_current, h_mock, p, zero(T_val))
    
    # A karakterisztikus mátrix-egyenlet: λ*A - (du - du_eq) = 0
    res = λ .* A_duals .- du
    
    # --- 2. JACOBI MÁTRIX KINYERÉSE (0 ALLOKÁCIÓS SMatrix) ---
    get_p(x, j) = ForwardDiff.partials(x, j)
    
    M_lambda = SMatrix{N, N, Complex{T_val}, N*N}(
        ntuple(k -> begin
            i = (k - 1) % N + 1
            j = (k - 1) ÷ N + 1
            r = res[i]
            complex(get_p(real(r), j), get_p(imag(r), j))
        end, Val(N*N))
    )
    
    # --- 3. FIXPONT ELLENŐRZÉS ---
    if verbosity > 0
        get_v(x) = ForwardDiff.value(x)
        # Itt eltávolítjuk a belső Tag-et, hogy megnézzük a maradékot
        val_eq = map(r -> complex(get_v(real(r)), get_v(imag(r))), res)
        
        # Rekurzív Float64 kinyerő segédfüggvény (ha λ kívülről duális volt)
        _raw_float(x::Float64) = x
        _raw_float(x::ForwardDiff.Dual) = _raw_float(ForwardDiff.value(x))
        
        sum_err = sum(abs2, val_eq)
        if _raw_float(sum_err) > 1e-6
            @warn "A megadott u_eq pontban a deriváltak nem nullák! Maradék (du): $(val_eq)"
        end
    end
    
    # Az SMatrix det() hívása loop-unrolled, explicit képlettel fut le!
    return det(M_lambda)
end


# ==============================================================================
# 3. TESZTELÉS ÉS BIZONYÍTÁS
# ==============================================================================
println("--- Feketedoboz Karakterisztikus Egyenlet Teszt ---")

p_test = (0.5, -0.1)
ω_test = 2.5
σ_test = -0.05
λ_test = σ_test + 1im * ω_test

# N_dof-ot átalakítjuk Val(2)-re a fordító optimalizációja miatt
D_chareq_wrapper(p1, p2, λ) = get_D_from_model(turning_model, λ, (p1,p2), Val(2), verbosity=0)
D_numerical = D_chareq_wrapper(p_test[1], p_test[2], λ_test)

D_analytical = D_chareq_analytical(p_test[1], p_test[2], λ_test)
println("Teszt pont: λ = $λ_test")
println("1. Analitikus D(λ) eredmény:  ", D_analytical)
println("2. Numerikus  D(λ) eredmény:  ", D_numerical)

hiba = abs(D_analytical - D_numerical)
println("\nA két módszer közötti különbség: ", hiba)

if hiba < 1e-10
    println("-> SIKER: A numerikus mátrix-kifejtés tökéletesen megegyezik az analitikussal!")
else
    println("-> HIBA: Eltérés a számításban.")
end

# ==============================================================================
# 4. TELJESÍTMÉNYMÉRÉS (BENCHMARK)
# ==============================================================================
using BenchmarkTools

println("\n--- Teljesítménymérés (Benchmark) ---")

println("\n1. Analitikus D(λ) benchmark:")
bench_analytical = @benchmark D_chareq_analytical($(p_test[1]), $(p_test[2]), $λ_test)
display(bench_analytical)

println("\n2. Numerikus (Feketedoboz) D(λ) benchmark:")
bench_numerical = @benchmark D_chareq_wrapper($(p_test[1]), $(p_test[2]), $λ_test)
display(bench_numerical)