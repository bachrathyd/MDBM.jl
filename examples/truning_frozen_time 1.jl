5 + 5
using MDBM

using Plots
using LaTeXStrings

# --- 0. System Definition ---

# Global constant for damping ratio
const ζ = 0.05

"""
Core characteristic equation.
Takes all inputs as Float64 and returns a stable Tuple{Float64, Float64}.
This is the single source of truth for the system's physics.
"""
function characteristic_eq_tuple(Ω::Float64, w::Float64, σ::Float64, ω_c::Float64)
    # Construct complex lambda (root)
    λ = complex(σ, ω_c)

    # Time delay
    τ = 2.0 * π / Ω

    # Characteristic Equation: D = λ² + 2ζλ + 1 + w(1 - e^(-λτ))
    D = λ^2 + 2.0 * ζ * λ + 1.0 + w * (1.0 - exp(-λ * τ))

    # Return a tuple of (Real, Imag)
    return (real(D), imag(D))
end


# --- 1. Stability Chart (σ = 0) ---

"""
Wrapper function for MDBM to find the stability boundary.
Fixes sigma = 0.0.
"""
function stability_eq(Ω::Float64, w::Float64, ω_c::Float64)
    return characteristic_eq_tuple(Ω, w, 0.0, ω_c)
end

println("--- 1. Computing Stability Chart (σ=0) ---")

# Define coarse axes (approx 10 points)
ax_Ω_stab = MDBM.Axis(range(0.1, 2.0, length=10), "Ω")
ax_w_stab = MDBM.Axis(range(0.0, 2.0, length=10), "w")
ax_ωc_stab = MDBM.Axis(range(0.5, 3.5, length=10), "ω_c")


# Set up the MDBM problem
prob_stab = MDBM_Problem(
    stability_eq,
    [ax_Ω_stab, ax_w_stab, ax_ωc_stab]
)

# Solve with higher iteration depth for refinement
solve!(prob_stab, 5)
sol_stab = getinterpolatedsolution(prob_stab)

# Plot the stability lobes
p1 = scatter(sol_stab[1], sol_stab[2],
    markersize=1.5, legend=false, color=:blue,
    xlabel=L"\Omega \text{ (Spindle Speed)}",
    ylabel=L"w \text{ (Depth of Cut)}",
    title="Stability Chart (Type-Stable)"
)
display(p1)


## --- 2. Root Plane Analysis (Fixed Ω, w) ---

# Fixed parameters for this section
const Ω_fix = 1.0
const w_fix = 0.5 # Chosen as an unstable point

"""
Wrapper for MDBM: Finds where Real(D) = 0.
Returns a 1-element tuple as required by MDBM.
"""
function root_eq_real(σ::Float64, ω_c::Float64)
    # We only return the first element (Real part)
    return (characteristic_eq_tuple(Ω_fix, w_fix, σ, ω_c)[1],)
end

"""
Wrapper for MDBM: Finds where Imag(D) = 0.
Returns a 1-element tuple.
"""
function root_eq_imag(σ::Float64, ω_c::Float64)
    # We only return the second element (Imag part)
    return (characteristic_eq_tuple(Ω_fix, w_fix, σ, ω_c)[2],)
end

"""
Wrapper for MDBM: Finds the roots (Real=0 AND Imag=0).
Returns the full 2-element tuple.
"""
function root_eq_both(σ::Float64, ω_c::Float64)
    return characteristic_eq_tuple(Ω_fix, w_fix, σ, ω_c)
end

println("--- 2. Analyzing Roots in the Complex Plane ---")

# Define coarse axes for the complex plane
ax_σ_plane = Axis(range(-0.5, 0.5, length=10), "σ")
ax_ωc_plane = Axis(range(0.0, 2.5, length=10), "ω_c")

# Problem A: Find Real(D) = 0 curves
prob_re = MDBM_Problem(root_eq_real, [ax_σ_plane, ax_ωc_plane])
solve!(prob_re, 5)
sol_re = getinterpolatedsolution(prob_re)

# Problem B: Find Imag(D) = 0 curves
prob_im = MDBM_Problem(root_eq_imag, [ax_σ_plane, ax_ωc_plane])
solve!(prob_im, 5)
sol_im = getinterpolatedsolution(prob_im)

# Problem C: Find Roots (Intersections)
prob_roots = MDBM_Problem(root_eq_both, [ax_σ_plane, ax_ωc_plane])
solve!(prob_roots, 5)
sol_roots = getinterpolatedsolution(prob_roots)

# Plot the results
p2 = plot(title="Root Locations for Ω=$(Ω_fix), w=$(w_fix)", xlabel=L"\sigma", ylabel=L"\omega_c")
scatter!(p2, sol_re[1], sol_re[2], markersize=1, color=:red, label="Real(D)=0")
scatter!(p2, sol_im[1], sol_im[2], markersize=1, color=:green, label="Imag(D)=0")
scatter!(p2, sol_roots[1], sol_roots[2], markersize=6, color=:black, shape=:star5, label="Roots")
display(p2)


# --- 3. Root Movement (Varying Ω) ---

"""
Wrapper for MDBM to find the root locus (movement).
w is fixed (using w_fix), Ω is a variable.
"""
function locus_eq(σ::Float64, ω_c::Float64, Ω::Float64)
    # w_fix is the constant depth of cut
    return characteristic_eq_tuple(Ω, w_fix, σ, ω_c)
end

println("--- 3. Root Movement (Varying Omega) ---")

# Define coarse axes for the 3D (σ, ω_c, Ω) search space
ax_σ_move = Axis(range(-0.5, 0.5, length=10), "σ")
ax_ωc_move = Axis(range(0.5, 2.0, length=10), "ω_c")
ax_Ω_move = Axis(range(0.2, 0.8, length=10), "Ω")

# Set up the MDBM problem
prob_move = MDBM_Problem(
    locus_eq,
    [ax_σ_move, ax_ωc_move, ax_Ω_move]
)

solve!(prob_move, 5)
sol_move = getinterpolatedsolution(prob_move)

# Plot 3D Locus
p3 = scatter(sol_move[1], sol_move[2], sol_move[3],
    zcolor=sol_move[3], # Color by Omega
    xlabel=L"\sigma", ylabel=L"\omega_c", zlabel=L"\Omega",
    title="Root Movement (Ω: 0.3 → 0.5)",
    markersize=2, camera=(40, 30), label=false,
    colorbar_title=L"\Omega"
)
display(p3)

# Plot 2D Projection
p4 = scatter(sol_move[1], sol_move[2],
    zcolor=sol_move[3],
    xlabel=L"\sigma", ylabel=L"\omega_c",
    title="Root Locus Projection (Color = Ω)",
    label=false, colorbar_title=L"\Omega"
)
display(p4)

println("--- All computations complete. ---")

 DT1 = connect(prob_move)
 edge2plot_xyz = [reduce(hcat, [i_sol[getindex.(DT1, 1)], i_sol[getindex.(DT1, 2)], fill(NaN, length(DT1))])'[:] for i_sol in sol_move]
 p5=plot(edge2plot_xyz..., linewidth=5, label="midpoints solution connected")
 
display(p5)




 function connected_edge_components(edges::Vector{Tuple{Int,Int}})
    # undirected adjacency
    adj = Dict{Int, Vector{Int}}()
    for (u, v) in edges
        push!(get!(adj, u, Int[]), v)
        push!(get!(adj, v, Int[]), u)
    end

    visited = Set{Int}()
    comps   = Vector{Vector{Tuple{Int,Int}}}()

    for start in keys(adj)
        start in visited && continue

        # --- BFS/DFS over nodes in this component
        stack = [start]
        push!(visited, start)
        nodes = Int[start]

        while !isempty(stack)
            u = pop!(stack)
            for v in adj[u]
                if !(v in visited)
                    push!(visited, v)
                    push!(stack, v)
                    push!(nodes, v)
                end
            end
        end

        # collect edges whose endpoints are in this node set
        nodeset = Set(nodes)
        comp_edges = [e for e in edges if (e[1] in nodeset && e[2] in nodeset)]
        push!(comps, comp_edges)
    end

    return comps
end
















comps = connected_edge_components(DT1)
length(comps)           # → should be 3 for your case

p8=plot()
for DTloc in comps
 edge2plot_xyz = [reduce(hcat, [i_sol[getindex.(DTloc, 1)], i_sol[getindex.(DTloc, 2)], fill(NaN, length(DTloc))])'[:] for i_sol in sol_move]
 plot!(edge2plot_xyz..., linewidth=5)
end

display(p8)


p9=plot()
for DTloc in comps
 edge2plot_xyz = [reduce(hcat, [i_sol[getindex.(DTloc, 1)], i_sol[getindex.(DTloc, 2)], fill(NaN, length(DTloc))])'[:] for i_sol in sol_move]
 plot!(edge2plot_xyz[1] ,edge2plot_xyz[2], linewidth=5)#./edge2plot_xyz[3].^0.6
end
plot!()

display(p9)

 #for DTloc=comps[1:5]

    DTloc=comps[5]
x,y,z = [i_sol[getindex.(DTloc, 1)] for i_sol in sol_move]
p   = sortperm(z)
zs  = z[p]
xs  = x[p]
ys  = y[p]

#edge2plot_xyz = [reduce(hcat, [i_sol[getindex.(DTloc, 1)], i_sol[getindex.(DTloc, 2)], fill(NaN, length(DTloc))])'[:] for i_sol in sol_move]
# plot(edge2plot_xyz[[3,1]]..., linewidth=5)
# plot!(edge2plot_xyz[[3,2]]..., linewidth=5)
p10=plot(zs,xs, linewidth=5)
plot!(zs,ys, linewidth=5)

#end

using Polynomials

deg = 3  # or 3–5, don’t go crazy

indrange=(zs .> 0.3) .& (zs .<=0.4)
px = fit(zs[indrange], xs[indrange], deg)   # x(ω)
py = fit(zs[indrange], ys[indrange], deg)   # y(ω)

#--------------- plotting the results ------------------

N      = 200
zmin   = minimum(zs) - 0.05  # extend if you want
zmax   = maximum(zs) + 0.05
zgrid  = range(zmin, zmax; length = N)

# using polynomial fit:
xgrid = px.(zgrid)
ygrid = py.(zgrid)

plot!(zgrid,xgrid, linewidth=2 , linestyle = :dash)
plot!(zgrid,ygrid, linewidth=2 , linestyle = :dash)

display(p10)


t=0:0.01:2*pi
zgrid=0.45 .+ 0.2.*sin.(t)
#zgrid=Omega0 (1 .+ RVA .*sin.(RVF*Omega0*t)
ygrid=px.(zgrid)
p11=plot(t,ygrid, linewidth=2 , linestyle = :dash)

display(p11)

# 
# using Interpolations
# 
# itp_x = interpolate((zs,), xs, Gridded(Linear()))
# itp_y = interpolate((zs,), ys, Gridded(Linear()))
# 
# # choose extrapolation behaviour for outside-range ω:
# ext_x = extrapolate(itp_x, Line())  # linear extrapolation
# ext_y = extrapolate(itp_y, Line())
# 
# 