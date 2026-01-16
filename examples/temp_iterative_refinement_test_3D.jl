5 + 5
using Revise
using MDBM
using GLMakie

using LinearAlgebra
GLMakie.closeall()
GLMakie.activate!(; title="2 parameters, codimension 1")
#-----------------------------

function foo_par2_codim1(x, y, z)
    pow = 7.0
    pow = 1.3
    #pow=0.75
    abs(x)^pow + abs(y)^pow + abs(z)^pow - 2.0^pow
    #sin(x*5)-y+z+ 0.1 * x * y
end
foo_par2_codim1(0.0, 0.0, 0.5)

#mymdbm = MDBM_Problem(foo_par2_codim1, [-0.01:3.0, -0.02:3.0, -0.02:3.0])
mymdbm = MDBM_Problem(foo_par2_codim1, [-1.06:0.5:3.0, -1.09:0.5:3.0, -1.1:0.5:3.0])
@time solve!(mymdbm, 2, interpolationorder=1, normp=90.0, ncubetolerance=-0.05)#number of refinements - increase it slightly to see smoother results 

##
nc_list = 1:size(mymdbm.ncubes, 1)


errv_s = [MDBM.getscaled_local_point(MDBM.ncube_error_vector(nc), nc, mymdbm.axes) for nc in mymdbm.ncubes]
err_norm = norm.(errv_s)
#nc_list = nc_list[err_norm.>sort(err_norm)[length(err_norm)÷10*5]]
nc_list = nc_list[err_norm.>=(minimum(err_norm)+maximum(err_norm))*0.5]#0.61803398875
#nc_list= nc_list[err_norm .> 0.01]


#axesextend!(mymdbm,1, -2.0:0.1:2)#both prepend and append automatically (overlapping values are eliminated)


# TODO: this is a problem, doubling must be done for only the cubes which hase size 1, and at the location where it is size 1 - this way i will not be able to tell the size (diffference) of the neighbouring n-cubes
# We doulging only in the directions where refinement is needed
nc_size_minimum_in_the_list = minimum([minimum(nc.size) for nc in mymdbm.ncubes[nc_list]])
println("-----------------------------------Mins size  $nc_size_minimum_in_the_list-------------------------------")
if nc_size_minimum_in_the_list < 2
    println("-----------------------------------Doubling...-------------------------------")
    MDBM.doubling!(mymdbm, [1, 2, 3])
end

@show length(nc_list)
MDBM.refinencubes!(mymdbm.ncubes, nc_list, [1, 2, 3])
MDBM.interpolate!(mymdbm, interpolationorder=1)

#@time checkneighbour!(mymdbm,verbosity=3);
#@time solve!(mymdbm, 1,interpolationorder=1)

# --- create once ---
if !@isdefined(fig)
    fig = Figure(size=(2300, 1350))

    # grid: left spans two rows; right has two stacked axes
    ax = GLMakie.Axis3(fig[1:2, 1], title="Result")
    ax_top = GLMakie.Axis(fig[1, 2], title="Error histogram", xlabel="‖error‖", ylabel="count")
    ax_bot = GLMakie.Axis(fig[2, 2], title="nc_size histogram", xlabel="nc_size", ylabel="count")

    display(fig)

    __hist_colors__ = Makie.wong_colors()
    __hist_colors__ = collect(cgrad([:blue, :red], 12, categorical=true))
    __run_idx__ = Ref(0)
end





MDBM.getcornerval(mymdbm.ncubes[1], mymdbm)
allcorners = MDBM.corner(mymdbm.ncubes, mymdbm.T01)

x = reduce(vcat, [[[mymdbm.axes[1][xy[1]] for xy in getindex(nc, [1, 2, 4, 3, 7, 8, 6, 5, 1])]..., NaN] for nc in allcorners])
y = reduce(vcat, [[[mymdbm.axes[2][xy[2]] for xy in getindex(nc, [1, 2, 4, 3, 7, 8, 6, 5, 1])]..., NaN] for nc in allcorners])
z = reduce(vcat, [[[mymdbm.axes[3][xy[3]] for xy in getindex(nc, [1, 2, 4, 3, 7, 8, 6, 5, 1])]..., NaN] for nc in allcorners])






# --- data prep (your variables assumed present) ---
errv_s = MDBM.ncube_error_vector.(mymdbm.ncubes)

errv_s = MDBM.getscaled_local_point.(errv_s, mymdbm.ncubes, Ref(mymdbm.axes))
err_norm = norm.(errv_s)
xy_sol = getinterpolatedsolution(mymdbm)           # (xs, ys)
ms = (err_norm ./ maximum(err_norm)) .^ 0.2 .* 10
nc_size = log2.([nc.size[1] for nc in mymdbm.ncubes])

# --- left: overwrite each run ---
empty!(ax)
#lines!(ax, x, y,z)
scatter!(ax, xy_sol...; color=err_norm, colormap=:reds)#, markersize=ms
#scatter!(ax, xy_sol...; color=nc_size, colormap=:reds)#, markersize=ms,

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

@show mymdbm







if false#true# false
    #--------------------------- Sub-cube interpolation----------------

    #calcuatin the sub-cubes interpolations stored in the mymdbm.ncubes[i].posinterp
    @time interpsubcubesolution!(mymdbm)#,normp=30,ncubetolerance=0.1)
    #extracting the resutls to from the 
    @time path2points = extract_paths(mymdbm)

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
##
fig