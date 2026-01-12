5 + 5
using Revise
using MDBM
using GLMakie

using LinearAlgebra
GLMakie.closeall()
GLMakie.activate!(; title="2 parameters, codimension 1")
#-----------------------------

function foo_par2_codim1(x, y)
    pow = 120.0
    pow = 20.0
    pow = 0.2
    abs(x)^pow + abs(y)^pow - 2.0^pow + 0.0 * x * y

    abs(x)^0.4 + abs((y))^1.4 - 1.0
    y - abs(sin(1 / x) / x)
    ((x^2.0 + y) - 1.0^2.0) * x #TODO: test with this function, too
end


function f_unwrapped(x, y)
    a, b = 0.2, 0.15
    r = sqrt(x^2 + y^2)
    th0 = atan2(y, x)
    k = round(Int, (r - a) / (2π * b) - th0 / (2π))
    th = th0 + 2π * k
    return r - (a + b * th)
end
mymdbm = MDBM_Problem(foo_par2_codim1, [-3.01:3.0, -3.02:3.0])
#mymdbm = MDBM_Problem(f_unwrapped, [-5.01:0.4:5.0, -5.01:0.4:5.0])
@time solve!(mymdbm, 3, interpolationorder=1, ncubetolerance=ncubetolerance=0.0)#number of refinements - increase it slightly to see smoother results 

# mymdbm = MDBM_Problem(foo_par2_codim1, [[-0.0,2.0], [-0.0,2.0]])
# @time solve!(mymdbm, 1,interpolationorder=1)#number of refinements - increase it slightly to see smoother results 

##
for _ in 1:3  
     nc_list = 1:size(mymdbm.ncubes, 1)


    errv_s = [MDBM.getscaled_local_point(MDBM.ncube_error_vector(nc), nc, mymdbm.axes) for nc in mymdbm.ncubes]
    err_norm = norm.(errv_s)
   # nc_list = nc_list[err_norm.>=sort(err_norm)[length(err_norm)÷10*7+1]]
    nc_list = nc_list[err_norm.>=(minimum(err_norm)+maximum(err_norm))/2]

   #    nc_list= nc_list[err_norm .> 0.0001]


    # TODO: this is a problem, doubling must be done for only he cubes which hase size 1, and at the location where it is size 1 - this way i will not be able to tell the size (diffference) of the neighbouring n-cubes
     nc_size_minimum_in_the_list = minimum([minimum(nc.size) for nc in mymdbm.ncubes[nc_list]])
     if nc_size_minimum_in_the_list<2
            MDBM.doubling!(mymdbm, [1, 2])
     end
     
    if length(nc_list) == 0
        break
    end
    @show length(nc_list)
    #MDBM.refinencubes!(mymdbm.ncubes, 1:size(mymdbm.ncubes, 1) ÷ 2, [1, 2])
    MDBM.refinencubes!(mymdbm.ncubes, nc_list, [1, 2])
    #MDBM.refinencubes!(mymdbm.ncubes, 1:size(mymdbm.ncubes, 1) , [1])
    #MDBM.refinencubes!(mymdbm.ncubes, 1:size(mymdbm.ncubes, 1) , [2])
    MDBM.interpolate!(mymdbm, interpolationorder=1,normp=10.0, ncubetolerance=ncubetolerance=0.4)
    #end
    #@time solve!(mymdbm, 1,interpolationorder=1)

    # --- create once ---
    if !@isdefined(fig)
        fig = Figure(size=(2300, 1350))

        # grid: left spans two rows; right has two stacked axes
        ax = GLMakie.Axis(fig[1:2, 1], title="Result")
        ax_top = GLMakie.Axis(fig[1, 2], title="Error histogram", xlabel="‖error‖", ylabel="count")
        ax_bot = GLMakie.Axis(fig[2, 2], title="nc_size histogram", xlabel="nc_size", ylabel="count")

        display(fig)

        __hist_colors__ = Makie.wong_colors()
        __run_idx__ = Ref(0)
    end





    MDBM.getcornerval(mymdbm.ncubes[1], mymdbm)
    allcorners = MDBM.corner(mymdbm.ncubes, mymdbm.T01)

    x = reduce(vcat, [[[mymdbm.axes[1][xy[1]] for xy in getindex(nc, [1, 2, 4, 3, 1])]..., NaN] for nc in allcorners])
    y = reduce(vcat, [[[mymdbm.axes[2][xy[2]] for xy in getindex(nc, [1, 2, 4, 3, 1])]..., NaN] for nc in allcorners])






    # --- data prep (your variables assumed present) ---
    errv_s = MDBM.ncube_error_vector.(mymdbm.ncubes)
    errv_s = MDBM.getscaled_local_point.(errv_s, mymdbm.ncubes, Ref(mymdbm.axes))
    err_norm = norm.(errv_s)
    xy_sol = getinterpolatedsolution(mymdbm)           # (xs, ys)
    ms = (err_norm ./ maximum(err_norm)) .^ 0.2 .* 10
    nc_size = log2.([nc.size[1] for nc in mymdbm.ncubes])

    # --- left: overwrite each run ---
    empty!(ax)
    lines!(ax, x, y)
    scatter!(ax, xy_sol...; color=err_norm, markersize=ms, colormap=:viridis)
    #scatter!(ax, xy_sol...; color=nc_size, markersize=ms, colormap=:viridis)

    # --- right-top: accumulate histograms across runs ---
    __run_idx__[] += 1
    c = __hist_colors__[mod1(__run_idx__[], length(__hist_colors__))]
    empty!(ax_top)
    hist!(ax_top, err_norm; bins=30, color=(c, 0.85), strokecolor=:black)

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




#     ##
#     #--------------------------- Sub-cube interpolation----------------
# 
#     #calcuatin the sub-cubes interpolations stored in the mymdbm.ncubes[i].posinterp
#     @time interpsubcubesolution!(mymdbm,normp=10,ncubetolerance=0.1)
#     #extracting the resutls to from the 
#     path2points = extract_paths(mymdbm)
# 
#     #extracting the unique points and plotting
#     puniq = unique(collect(Iterators.flatten(Iterators.flatten(path2points))))
#     scatter!(ax, getindex.(puniq, 1), getindex.(puniq, 2), label="subface - solution")
# 
# 
# 
#     #exctracing the simplexes for each ncube
#     flatened_path2points = collect(Iterators.flatten(path2points))
#     #eliminating the points with less than 2 points (caused by fininte precision)
#     true_truflatened_path2points = flatened_path2points[length.(flatened_path2points).==2]
#     #plotting the lines between the points
#     lines2plot = [(Point2f(ploc[1]), Point2f(ploc[2])) for ploc in true_truflatened_path2points]
#     #empty!(ax)
#     linesegments!(ax, lines2plot, label="subface - connection")
# 
# 
# #    display(GLMakie.Screen(), fig)




end
##



checkneighbour!(mymdbm)

@show mymdbm


