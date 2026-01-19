5 + 5
using Revise
using MDBM
using GLMakie

using LinearAlgebra
GLMakie.closeall()

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

@time solve!(turning_mdbm_prob, 4, verbosity=0, checkneighbourNum=2)# check neighbour at every iteration
##
@time for _ in 1:3
nc_list = 1:size(turning_mdbm_prob.ncubes, 1)


errv_s = [MDBM.getscaled_local_point(MDBM.ncube_error_vector(nc), nc, turning_mdbm_prob.axes) for nc in turning_mdbm_prob.ncubes]
err_norm = norm.(errv_s)
nc_list = nc_list[err_norm.>sort(err_norm)[length(err_norm)÷10*5]]
#nc_list = nc_list[err_norm.>=(minimum(err_norm)+maximum(err_norm))/2]
#nc_list= nc_list[err_norm .> 0.01]
@show length(nc_list)

# TODO: this is a problem, doubling must be done for only he cubes which hase size 1, and at the location where it is size 1 - this way i will not be able to tell the size (diffference) of the neighbouring n-cubes
nc_size_minimum_in_the_list = minimum([minimum(nc.size) for nc in turning_mdbm_prob.ncubes[nc_list]])
println("-----------------------------------Mins size  $nc_size_minimum_in_the_list-------------------------------")
if nc_size_minimum_in_the_list < 2
    println("-----------------------------------Doubling...-------------------------------")
    MDBM.doubling!(turning_mdbm_prob, [1, 2, 3])
end


MDBM.refinencubes!(turning_mdbm_prob.ncubes, nc_list, [1, 2, 3])
MDBM.interpolate!(turning_mdbm_prob, interpolationorder=1)

#@time solve!(turning_mdbm_prob, 1,interpolationorder=1)

#checkneighbour!(turning_mdbm_prob,maxiteration=100)
#@show turning_mdbm_prob

# --- create once ---
if !@isdefined(fig)
    fig = Figure(size=(1800, 1050))

    # grid: left spans two rows; right has two stacked axes
    ax3D = GLMakie.Axis3(fig[1, 1], title="Result")
    ax2D = GLMakie.Axis(fig[2, 1], title="Result")
    ax_top = GLMakie.Axis(fig[1, 2], title="Error histogram", xlabel="‖error‖", ylabel="count")
    ax_bot = GLMakie.Axis(fig[2, 2], title="nc_size histogram", xlabel="nc_size", ylabel="count")

    display(fig)

    __hist_colors__ = Makie.wong_colors()
    __hist_colors__ = collect(cgrad([:blue, :red], 12, categorical=true))
    __run_idx__ = Ref(0)
end





MDBM.getcornerval(turning_mdbm_prob.ncubes[1], turning_mdbm_prob)
allcorners = MDBM.corner(turning_mdbm_prob.ncubes, turning_mdbm_prob.T01)

x = reduce(vcat, [[[turning_mdbm_prob.axes[1][xy[1]] for xy in getindex(nc, [1, 2, 4, 3, 7, 8, 6, 5, 1])]..., NaN] for nc in allcorners])
y = reduce(vcat, [[[turning_mdbm_prob.axes[2][xy[2]] for xy in getindex(nc, [1, 2, 4, 3, 7, 8, 6, 5, 1])]..., NaN] for nc in allcorners])
z = reduce(vcat, [[[turning_mdbm_prob.axes[3][xy[3]] for xy in getindex(nc, [1, 2, 4, 3, 7, 8, 6, 5, 1])]..., NaN] for nc in allcorners])






# --- data prep (your variables assumed present) ---
errv_s = MDBM.ncube_error_vector.(turning_mdbm_prob.ncubes)

errv_s = MDBM.getscaled_local_point.(errv_s, turning_mdbm_prob.ncubes, Ref(turning_mdbm_prob.axes))
err_norm = norm.(errv_s)
xy_sol = getinterpolatedsolution(turning_mdbm_prob)           # (xs, ys)
ms = (err_norm ./ maximum(err_norm)) .^ 0.2 .* 10
nc_size = log2.([nc.size[1] for nc in turning_mdbm_prob.ncubes])

# --- left: overwrite each run ---
empty!(ax3D)
#lines!(ax3D, x, y,z)
scatter!(ax3D, xy_sol...; color=err_norm, colormap=:viridis)#, markersize=ms
#scatter!(ax, xy_sol...; color=nc_size, colormap=:reds)#, markersize=ms,

empty!(ax2D)
lines!(ax2D, x, y)
#scatter!(ax2D, xy_sol[1], xy_sol[2]; color=nc_size, colormap=:reds)#, markersize=ms
# --- right-top: accumulate histograms across runs ---
__run_idx__[] += 1
c = __hist_colors__[mod1(__run_idx__[], length(__hist_colors__))]
empty!(ax_top)
hist!(ax_top, err_norm; bins=30, color=(c, 0.45), strokecolor=:black)

# --- right-bottom: overwrite each run with nc_size histogram ---
empty!(ax_bot)
# integer-centered bins for nicer bars if nc_size are integers
if !isempty(nc_size)
    lo, hi = minimum(nc_size), maximum(nc_size)
    edges = (lo-0.5):(hi+0.5)
    hist!(ax_bot, nc_size; bins=edges, color=(:gray, 0.6), strokecolor=:black)
end

fig

ax_top.yscale = log10
ax_bot.yscale = log10

#@show turning_mdbm_prob







if false #true#
    #--------------------------- Sub-cube interpolation----------------

    #calcuatin the sub-cubes interpolations stored in the turning_mdbm_prob.ncubes[i].posinterp
    @time interpsubcubesolution!(turning_mdbm_prob)#,normp=30,ncubetolerance=0.1)
    #extracting the resutls to from the 
    @time path2points = extract_paths(turning_mdbm_prob)

    #extracting the unique points and plotting
    puniq = unique(collect(Iterators.flatten(Iterators.flatten(path2points))))
    scatter!(ax, getindex.(puniq, 1), getindex.(puniq, 2), getindex.(puniq, 3), markersize=5, color=:green, label="subface - solution")



    #exctracing the simplexes for each ncube
    flatened_path2points = collect(Iterators.flatten(path2points))
    #eliminating the points with less than 2 points (caused by fininte precision)
    true_truflatened_path2points = flatened_path2points[length.(flatened_path2points).==3]


    #plotting the lines between the points
    n_faces = reshape(1:(3*length(true_truflatened_path2points)), (3, length(true_truflatened_path2points)))'
    vertices_mat = hcat(Iterators.flatten(true_truflatened_path2points)...)


    empty!(ax)
    mesh!(ax, vertices_mat, n_faces, alpha=0.995, label="subface - local simplex")

end
fig
end
##
println("Checking neighbours...end")
@time checkneighbour!(turning_mdbm_prob,maxiteration=15,verbosity=0);


nothing
#3.7  - with extra variable - first run
#4.5 - 5.3 sec - with extra variable - further runs

#5.2 - 3.5 sec - without extra variable - first run
#4.4 - 5.3 sec - without extra variable - further runs

# @profview  checkneighbour!(turning_mdbm_prob,maxiteration=20,verbosity=0);

# bulk - slow

#MDBM.getcornerval(turning_mdbm_prob.ncubes[1], turning_mdbm_prob)
#allcorners = MDBM.corner(turning_mdbm_prob.ncubes, turning_mdbm_prob.T01)
#
#x = reduce(vcat, [[[turning_mdbm_prob.axes[1][xy[1]] for xy in getindex(nc, [1, 2, 4, 3, 7, 8, 6, 5, 1])]..., NaN] for nc in allcorners])
#y = reduce(vcat, [[[turning_mdbm_prob.axes[2][xy[2]] for xy in getindex(nc, [1, 2, 4, 3, 7, 8, 6, 5, 1])]..., NaN] for nc in allcorners])
#z = reduce(vcat, [[[turning_mdbm_prob.axes[3][xy[3]] for xy in getindex(nc, [1, 2, 4, 3, 7, 8, 6, 5, 1])]..., NaN] for nc in allcorners])
#
#empty!(ax2D)
#lines!(ax2D, x, y)
#
#xy_sol = getinterpolatedsolution(turning_mdbm_prob)           # (xs, ys)
## --- left: overwrite each run ---
#empty!(ax3D)
##lines!(ax3D, x, y,z)
#scatter!(ax3D, xy_sol...)#, markersize=ms
##scatter!(ax, xy_sol...; color=nc_size, colormap=:reds)#, markersize=ms,
#
#fig
#