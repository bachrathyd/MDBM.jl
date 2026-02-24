
using MDBM
using GLMakie
using Statistics
using DSP # for unwrapping

using MAT

GLMakie.closeall()
GLMakie.activate!(; title="Encircle test")
cmapall = :rainbow

using MAT
using Interpolations

# 1. Load the .mat file
file = matopen(".\\examples\\depdin_FRFs.mat")
f_vector = read(file, "f")      # 1x11001 double
frf_matrix = read(file, "FRFs") # 29x11001 complex double
close(file)

om_scale = 20
stiff_scale = 22.0
Ndrop = 10
f_vector = f_vector[1, Ndrop:1000] ./ om_scale
first_frf = stiff_scale .* (frf_matrix[14, Ndrop:1000] ./ f_vector .^ 2)#Note, the sign change!
f_vector = vcat([-100000000, -f_vector[end] - 1], -f_vector[end:-1:1], f_vector, [f_vector[end] + 1, 100000000])
first_frf = vcat([0.0, 0.0], conj.(first_frf[end:-1:1]), first_frf, [0.0, 0.0])
k = 10
H_smooth = filtfilt(ones(k) / k, [1.0], first_frf)
const H_acc = LinearInterpolation(f_vector[:], H_smooth)

omv = LinRange(-1.012, 15.0, 70000)
#H_scaled(ω::Float64)::ComplexF64 = H_acc(ω)

f = Figure(size=(1700, 600))
ax1 = GLMakie.Axis(f[1, 1])
scatter!(omv, abs.(H_acc.(omv)))
scatter!(omv, real.(H_acc.(omv)))
scatter!(omv, imag.(H_acc.(omv)))
ζ = 0.02
λ = 1im * omv
D = 1.0 ./ (λ .^ 2 .+ 2 * ζ .* λ .+ 1.0)
scatter!(omv, abs.(D))
scatter!(omv, real.(D))
f


##
#  # val = H_acc(-150.5) 

# println("Interpolated value at 150.5 Hz: ", val)
# -----------------------------



#function D_chareq(p, d, ω)::Tuple{Float64,Float64}
#    λ = 1im * ω
#    τ = 0.5
#    ζ = 0.02
#    D = (0.03 * λ^4 + λ^2 + 2 * ζ * λ + 1 + p * exp(-τ * λ) + d * λ * exp(-τ * λ)) #/ (abs(λ)^2+ 1)
#    return real(D), imag(D)
#end
#n_power_max = 4

# function D_chareq(p::Float64, d::Float64, ω::Float64)::Tuple{Float64,Float64}
#     λ = 1im * ω
#     τ = 0.5
#     ζ = 0.02
#     #H=1/(λ^2 + 2 * ζ * λ + 1)
#     H = H_acc(ω)
#     D = (1 / H  + (p * exp(-τ * λ) + d * λ * exp(-τ * λ))) #/ (abs(λ)^2+ 1)
#     return real(D), imag(D)
# end
# n_power_max = 2


# function D_chareq(δ, b, ω)::Tuple{Float64,Float64}
#     λ = 1im * ω
#     τ = 2π
#     ζ = 0.02
#     D = (λ^2 + 2 * ζ * λ + δ + b * exp(-τ * λ)) / (abs(λ)^2 + 1)
#     return real(D), imag(D)
# end
# n_power_max = 2

function D_chareq(invΩ, w, ω)::Tuple{Float64,Float64}
    λ = 1im * ω
    τ = 2π*invΩ
    ζ = 0.01
    
    #H=1/(λ^2 + 2 * ζ * λ + 1)
    H = H_acc(ω)
    D = (1 / H + 1 + w*(1.0- exp(-τ * λ))) / (abs(λ)^2 + 1)
    return real(D), imag(D)
end
n_power_max = 2

#D_chareq(1.0, 1.0, 1.0)
Pv = LinRange(-2.01, 4.0, 20)
#Pv = LinRange(0.21, 4.0, 20)
Dv = LinRange(-2.01, 5.0, 20)
omv = [-100, LinRange(-3.0, 15.0, 40)..., 100, 1000]
#omv = [LinRange(0, 15.0, 70)..., (LinRange(15.1, 100, 30) .^ 2)...]
omv = LinRange(-0.012, 5.0, 70)
omv = [LinRange(-0.012, 5.0, 70)..., (LinRange(5.01, 10, 10) .^ 2)...]

omextra = [0.0]##collect(LinRange(0, 5, 35))# [0.0, 100.0];
omextra = sort(vcat([0.0], omv))

mymdbm = MDBM_Problem(D_chareq, [Pv, Dv, omv])


#@time solve!(mymdbm, 3, verbosity=1) #number of refinements - increase it slightly to see smoother results 
Niter = 5
@time MDBM.solve!(mymdbm, Niter, verbosity=0, checkneighbourNum=1, doThreadprecomp=true)

#@show mymdbm

f = Figure(size=(1700, 600))
ax1 = GLMakie.Axis3(f[1, 1])
# # n-cube interpolation
xyz_sol = getinterpolatedsolution(mymdbm)
# scatter!(ax1, xyz_sol..., markersize=6, color=:red, marker='x', strokewidth=3, label="solution")

# connecting and plotting the "mindpoints" of the n-cubes
DT1 = MDBM.connect(mymdbm)
edge2plot_xyz = [reduce(hcat, [i_sol[getindex.(DT1, 1)], i_sol[getindex.(DT1, 2)], fill(NaN, length(DT1))])'[:] for i_sol in xyz_sol]
lines!(ax1, edge2plot_xyz..., linewidth=3, label="midpoints solution connected")
display(f)

## MDBM direct Encircle number
println(" ----------------- MDBM direct Encircle number -----------------")
const om_fix = sort(mymdbm.axes[3], rev=false)
function Encirct(P, D)::Float64
    D_full = D_chareq.(Ref(P), Ref(D), om_fix)
    D_full_Comp = [D[1] + 1im * D[2] for D in D_full]
    angles = angle.(D_full_Comp)
    unwrapped_angles = unwrap(angles)
    total_phase_change = unwrapped_angles[end] - unwrapped_angles[1]
    N_BF = -total_phase_change / (2π) + n_power_max / 2
    return N_BF - 0.5
end

mymdbm_2D = MDBM_Problem(Encirct, [Pv, Dv])

@time MDBM.solve!(mymdbm_2D, Niter, verbosity=0, checkneighbourNum=1, doThreadprecomp=true)

xyz_sol = getinterpolatedsolution(mymdbm_2D)
# scatter!(ax1, xyz_sol..., markersize=6, color=:red, marker='x', strokewidth=3, label="solution")
DT1 = MDBM.connect(mymdbm_2D)
edge2plot_xyz = [reduce(hcat, [i_sol[getindex.(DT1, 1)], i_sol[getindex.(DT1, 2)], fill(NaN, length(DT1))])'[:] for i_sol in xyz_sol]
lines!(ax1, edge2plot_xyz..., linewidth=3, label="midpoints solution connected")
display(f)
## ----------------- BF ----------------
println(" ----------------- Start BF -----------------")
Pv = mymdbm.axes[1]#LinRange(-2.01, 4.0, 300)
Dv = mymdbm.axes[2]#LinRange(-2.01, 5.0, 300)

#omv=collect(omv)
omv = mymdbm.axes[3]#[-50, LinRange(-15.0, 15.0, 140)..., 50, 1000, 10000]
sort!(omv, rev=false)
N_BF = zeros(length(Pv), length(Dv))
sum_re = zeros(length(Pv), length(Dv))
@time for (pi, P) in enumerate(Pv)
    for (di, D) in enumerate(Dv)
        D_full = D_chareq.(Ref(P), Ref(D), omv)
        D_full_Comp = [D[1] + 1im * D[2] for D in D_full]
        angles = angle.(D_full_Comp)
        unwrapped_angles = unwrap(angles)
        total_phase_change = unwrapped_angles[end] - unwrapped_angles[1]

        N_BF[pi, di] = -total_phase_change / (2π) + n_power_max / 2
    end
end

#delete!(ax_2D)
ax_2D = GLMakie.Axis(f[1, 2])
sf = surface!(ax_2D, Pv, Dv, N_BF, colormap=cmapall, colorrange=(0, maximum(N_BF)), shading=false)
Colorbar(f[1, 2][1, 2], sf)#, vertical=false)

#contour!(ax_2D, Pv, Dv, N_BF, levels=[0.5], linewidth=4)
display(f)








## ------------MDBM------------------------
println(" ----------------- Start MDBM -----------------")



D_re_im = getevaluatedfunctionvalues(mymdbm)
D_comp = [D[1] + 1im * D[2] for D in D_re_im]
xyz_val = getevaluatedpoints(mymdbm)




xy_points = [[x, y] for (x, y, z) in zip(xyz_val...,)]
xy_p_uniq = unique(xy_points)


Ncirc = zeros(length(xy_p_uniq))
# using LinearAlgebra
# i = argmin(norm.(xy_p_uniq .- Ref([0.192, 0.189])))
P = 0.35
D = 0.18
# i = argmin(norm.(xy_p_uniq .- Ref([P, D])))
# xy_loc = xy_p_uniq[i]


# #.............BF...................
## ax_4 = GLMakie.Axis(f[1, 4])
# omv = LinRange(-150.0, 150.0,5140)
# D_full = D_chareq.(Ref(P), Ref(D), omv)
# D_full_Comp = [D[1] + 1im * D[2] for D in D_full]
# scatter!(ax_4, real.(D_full_Comp), imag(D_full_Comp))

# #.............MDBM...................
@time for (i, xy_loc) in enumerate(xy_p_uniq)
    # if mod(i, 100) == 0
    #     println(i / length(xy_p_uniq))
    # end
    #xy_loc=xy_p_uniq[597]
    indexing = Ref(xy_loc) .== xy_points
    om_loc = xyz_val[3][indexing]
    Dloc = D_comp[indexing]
    #@show length(Dloc)

    om_loc = vcat(om_loc, omextra)
    D_extra = D_chareq.(Ref(xy_loc[1]), Ref(xy_loc[2]), omextra)
    Dloc = vcat(Dloc, [D[1] + 1im * D[2] for D in D_extra])

    # Assuming you have your data: om_loc and Dloc (Complex numbers)

    # 1. Create the symmetric data (Negative frequencies)
    # Since D(-ω) = conj(D(ω)) for real-coefficient systems
    # 1. Combine the original and its complex conjugate mirror
    om_merged = vcat(om_loc, -om_loc)
    D_merged = vcat(Dloc, conj.(Dloc))

    # 2. Get the sorting indices to keep om and D synchronized
    p = sortperm(om_merged, rev=true)

    # 3. Apply sorting and filter out duplicate omega values
    # 'unique' on the indices ensures we don't have two points at the same frequency
    mask = unique(i -> om_merged[i], p)

    om_full = om_merged[mask]
    D_full = D_merged[mask]


    # ax_4 = GLMakie.Axis(f[1, 4])
    # scatter!(ax_4, real.(D_full), imag.(D_full))

    # 2. Calculate the phase of each point
    angles = angle.(D_full)

    # 3. Calculate the change in angle between successive points
    # Unwrap handles the 2π jumps to ensure a continuous phase curve
    unwrapped_angles = unwrap(angles)

    # 4. Total encirclements N = (Total Phase Change) / (2π)
    total_phase_change = unwrapped_angles[end] - unwrapped_angles[1]
    N = total_phase_change / (2π)

    # println("Total Encirclements of Zero: ", round(N, digits=2))


    Ncirc[i] = N + n_power_max / 2


end



ax_2D_2 = GLMakie.Axis(f[1, 3])
sf = scatter!(ax_2D_2, getindex.(xy_p_uniq, 1), getindex.(xy_p_uniq, 2),
    color=Ncirc, colormap=cmapall, colorrange=(0, maximum(Ncirc)), markersize=10, label="evaluated")
Colorbar(f[1, 3][1, 2], sf)#, vertical=false)

xyz_sol = getinterpolatedsolution(mymdbm)
DT1 = MDBM.connect(mymdbm)
edge2plot_xyz = [reduce(hcat, [i_sol[getindex.(DT1, 1)], i_sol[getindex.(DT1, 2)], fill(NaN, length(DT1))])'[:] for i_sol in xyz_sol]
lines!(edge2plot_xyz..., linewidth=5, label="midpoints solution connected")
display(f)

#delete!(ax_3D_2)
ax_3D_2 = GLMakie.Axis(f[1, 4])
using DelaunayTriangulation
using GeometryBasics
using Statistics
x = getindex.(xy_p_uniq, 1)
y = getindex.(xy_p_uniq, 2)
z = Ncirc
# x = vcat(getindex.(xy_p_uniq, 1),xyz_sol[1])
# y =  vcat(getindex.(xy_p_uniq, 2),xyz_sol[2])
# z =  vcat(Ncirc,xyz_sol[3].*0.0)
#
# 1. Triangulate
@time tri = triangulate([x'; y'])

# 2. Convert triangles to a format Makie understands (Faces)
# We use each_solid_triangle to avoid any ghost triangles from the boundary
faces = [TriangleFace(t[1], t[2], t[3]) for t in each_solid_triangle(tri)]

# 3. Combine points (x, y, z) into 3D points
points = Point3f.(x, y, z)

# 4. Create the mesh object
m = GeometryBasics.mesh(points, faces)

# 5. Plot the mesh (Notice we don't pass 'indices' here, they are inside 'm')
mh = mesh!(ax_3D_2, m, color=z, colormap=cmapall, shading=false, colorrange=(0, maximum(Ncirc)))
Colorbar(f[1, 4][1, 2], mh)#, vertical=false)

# 6. Optional: Add the wireframe to see the triangulation
wireframe!(ax_3D_2, m, color=(:black, 0.1))
display(f)



#### ****************************************************************************

println(" ----------------- Start Triangulation -----------------")
using DelaunayTriangulation, GLMakie, Statistics, GeometryBasics
@time begin
    # 1. Setup points as Point2f
    pts_solution = [Point2f(xyz_sol[1][i], xyz_sol[2][i]) for i in 1:length(xyz_sol[1])]
    pts_sampled = [Point2f(p[1], p[2]) for p in xy_p_uniq]
    all_pts = vcat(pts_solution, pts_sampled)
    n_sol = length(pts_solution)

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
        global face_count
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
end
# 4. Plotting
ax_5 = GLMakie.Axis(f[1, 5], title="Stability Regions (Ncirc)")

# Construct the explicit mesh
m_obj = GeometryBasics.Mesh(mesh_points, mesh_faces)

# Plot with shading disabled to avoid the specular crash
mh = mesh!(ax_5, m_obj,
    color=round.(mesh_colors),#round.
    colormap=cmapall, colorrange=(0, maximum(Ncirc)),
    shading=NoShading) # Prevents the ComputePipeline error

Colorbar(f[1, 5][1, 2], mh)#, vertical=false)
# Overlay boundaries
lines!(ax_5, edge2plot_xyz..., color=:white, linewidth=1)


#scatter!(ax_5, xyz_sol, xyz_sol[2], markersize=10,  marker='.')
display(f)