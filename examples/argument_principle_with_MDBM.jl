#Argument Principle (Nyquist stability criterion)
# based on a soltuion of a characeristic equation 
using MDBM
using Statistics
using DSP # for unwrapping

using Interpolations


using DelaunayTriangulation, Statistics, GeometryBasics



"""
    calculate_encirclement_number(D_values, omega_values; n_power_max=0)

Calculates the number of encirclements of the origin for a given list of complex values `D_values`
corresponding to a list of frequencies `omega_values`. This is often used in Nyquist stability analysis.

# Arguments
- `D_values::AbstractVector{<:Complex}`: A vector of complex numbers representing the path in the complex plane.
- `omega_values::AbstractVector{<:Real}`: The corresponding frequencies for each point in `D_values`.
- `n_power_max::Integer=0`: A parameter related to the number of poles of the system at infinity,
  used to complete the stability criterion calculation.

# Returns
- `Float64`: The calculated number of encirclements, adjusted by `n_power_max`.
"""
function calculate_encirclement_number(D_values::AbstractVector{<:Complex}, omega_values::AbstractVector{<:Real}; n_power_max::Integer=0)
    if length(D_values) != length(omega_values)
        error("D_values and omega_values must have the same length.")
    end
    if length(D_values) < 2
        return Float64(n_power_max / 2) # Not enough points for phase change
    end

    # Sort by omega to trace the path correctly
    p = sortperm(omega_values)
    omega_sorted = omega_values[p]
    D_sorted = D_values[p]

    # Filter unique omega values, keeping the point for the last occurrence in the sorted list
    mask = unique(i -> omega_sorted[i], 1:length(omega_sorted))
    omega_unique = omega_sorted[mask]
    D_unique = D_sorted[mask]

    if length(D_unique) < 2
        return Float64(n_power_max / 2) # Not enough unique points
    end

    # Calculate phase and unwrap to get continuous phase change
    angles = angle.(D_unique)
    unwrapped_angles = unwrap(angles)
    total_phase_change = unwrapped_angles[end] - unwrapped_angles[1]

    # Number of encirclements (negative sign for standard control theory convention)
    N = -total_phase_change / (2π)

    # Add contribution from poles at infinity
    return N + n_power_max / 2
end

"""
    calculate_encirclement_number(D_chareq_handle, omega_values, params...; n_power_max=0)

Calculates the number of encirclements by evaluating a characteristic equation `D_chareq_handle`
over a range of frequencies `omega_values` with given parameters `params`.

# Arguments
- `D_chareq_handle::Function`: A function with signature `(params..., ω::Real)` that returns a complex number.
- `omega_values::AbstractVector{<:Real}`: The frequencies to evaluate `D_chareq_handle` at.
- `params...`: The parameters to pass to `D_chareq_handle`.
- `n_power_max::Integer=0`: Number of poles at infinity.

# Returns
- `Float64`: The calculated number of encirclements.
"""
function calculate_encirclement_number(D_chareq_handle::Function, omega_values::AbstractVector{<:Real}, params...; n_power_max::Integer=0)
    # Evaluate the characteristic equation at each frequency
    D_complex = [D_chareq_handle(params..., ω) for ω in omega_values]

    # Use the other method to do the actual calculation
    return calculate_encirclement_number(D_complex, omega_values; n_power_max=n_power_max)
end


function argument_principle_with_MDBM(D_chareq,mdbm, omega_coars, n_power_max)
    # 1. Get data from MDBM
    D_re_im = getevaluatedfunctionvalues(mdbm)
    D_comp = [D[1] + 1im * D[2] for D in D_re_im]

    xyz_val = getevaluatedpoints(mdbm)
    p1_vals, p2_vals, omega_vals = xyz_val

    # 2. Group results by parameter pair for efficiency
    grouped_results = Dict{Tuple{Float64,Float64},Vector{Tuple{Float64,ComplexF64}}}()
    sizehint!(grouped_results, length(p1_vals) ÷ 5) # Pre-allocate memory for the dictionary
    for i in 1:length(p1_vals)
        pair = (p1_vals[i], p2_vals[i])
        datapoint = (omega_vals[i], D_comp[i])
        if !haskey(grouped_results, pair)
            grouped_results[pair] = []
        end
        push!(grouped_results[pair], datapoint)
    end

    xy_p_uniq_tuples = keys(grouped_results) |> collect
    xy_p_uniq = [[p[1], p[2]] for p in xy_p_uniq_tuples] # Keep original format for plotting
    Ncirc = zeros(length(xy_p_uniq))

    # 3. Loop over unique parameter pairs and calculate Ncirc
    @time for (i, xy_loc_tuple) in enumerate(xy_p_uniq_tuples)
        # a. Get MDBM results for this pair
        mdbm_results = grouped_results[xy_loc_tuple]
        om_mdbm = getindex.(mdbm_results, 1)
        D_mdbm = getindex.(mdbm_results, 2)

        # b. Calculate D for the extra low-res points
        D_extra = D_chareq.(Ref(xy_loc_tuple[1]), Ref(xy_loc_tuple[2]), omega_coars)

        # c. Combine MDBM, extra, and symmetric points
        om_combined = vcat(om_mdbm, omega_coars)
        D_combined = vcat(D_mdbm, D_extra)
        om_full = vcat(om_combined, -om_combined)
        D_full = vcat(D_combined, conj.(D_combined))

        # d. Use the dedicated function to calculate encirclement number
        Ncirc[i] = calculate_encirclement_number(D_full, om_full; n_power_max=n_power_max)
    end

    return xy_p_uniq, Ncirc
end


function trinangulation_of_MDBM_results(mdbm, xy_p_uniq, Ncirc)

    xyz_sol = getinterpolatedsolution(mdbm)
    # 1. Setup points as Point2f
    pts_solution = [Point2f(xyz_sol[1][i], xyz_sol[2][i]) for i in 1:length(xyz_sol[1])]
    pts_sampled = [Point2f(p[1], p[2]) for p in xy_p_uniq]
    all_pts = vcat(pts_solution, pts_sampled)
    n_sol = length(pts_solution)

    DT1 = MDBM.connect(mdbm)
    edge2plot_xyz = [reduce(hcat, [i_sol[getindex.(DT1, 1)], i_sol[getindex.(DT1, 2)], fill(NaN, length(DT1))])'[:] for i_sol in xyz_sol]

    # 2. Triangulate with constraints
    const_edges = Set{Tuple{Int,Int}}(DT1)
    tri = triangulate(all_pts, segments=const_edges)

    # 3. Build a "Decomposed" Mesh (One color per triangle)
    # Instead of sharing vertices, we create 3 unique vertices per triangle
    # This ensures the color buffer length matches the vertex buffer length.
    mesh_points = Point2f[]
    mesh_faces = GeometryBasics.TriangleFace{Int}[]
    mesh_colors = Float32[]
    face_count = 0

    for T in each_triangle(tri)
        #global face_count
        u, v, w = T
        # Skip ghost triangles
        u > 0 && v > 0 && w > 0 || continue

        # Calculate average Ncirc for this triangle
        vals = [Ncirc[idx-n_sol] for idx in (u, v, w) if idx > n_sol]

        if !isempty(vals)
            avg_c = Float32(mean(vals))

            # Add the 3 vertices specifically for this face
            push!(mesh_points, all_pts[u], all_pts[v], all_pts[w])

            # Add the face using the new indices
            push!(mesh_faces, GeometryBasics.TriangleFace(face_count * 3 + 1, face_count * 3 + 2, face_count * 3 + 3))

            # Add the same color for all 3 vertices of this face
            append!(mesh_colors, [avg_c, avg_c, avg_c])

            face_count += 1
        end
    end
    return mesh_points, mesh_faces, mesh_colors,edge2plot_xyz
end

function argument_principle_solver_with_MDBM(Dchar_foo, axlist, n_power_max, ω_coars)
    @warn "Don't be lazy, adjust the settings of each steps for your needs!"
    function Dchar_foo_ReIm(p...)::Tuple{Float64,Float64}
        D = Dchar_foo(p...,)
        return real(D), imag(D)
    end
    # ----------------- Bifurcation curves with MDBM solution -----------------
    println(" ----------------- Bifurcation curves with MDBM solution -----------------")
    boundary_mdbm = MDBM_Problem(Dchar_foo_ReIm, axlist)

    Niter = 4
    @time MDBM.solve!(boundary_mdbm, Niter, verbosity=1, checkneighbourNum=2, doThreadprecomp=true, normp=10.0, ncubetolerance=0.6)

    ## ------------MDBM - updated Encircle number ------------------------
    println(" ----------------- MDBM - updated Encircle number -----------------")
    xy_p_uniq, Ncirc = argument_principle_with_MDBM(Dchar_foo,boundary_mdbm, ω_coars, n_power_max)

    ## step - 3 - triangulation for beautiful high res. visualization of the stability chart
    println(" ----------------- triangulation for beautiful high res. visualization of the stability chart -----------------")
    @time mesh_points, mesh_faces, mesh_colors,edge2plot_xyz = trinangulation_of_MDBM_results(boundary_mdbm, xy_p_uniq, Ncirc)

    return boundary_mdbm, xy_p_uniq, Ncirc, mesh_points, mesh_faces, mesh_colors,edge2plot_xyz
end