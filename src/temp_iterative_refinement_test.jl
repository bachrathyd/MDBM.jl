5 + 5
using Revise
using MDBM
using GLMakie

GLMakie.closeall()
GLMakie.activate!(; title="2 parameters, codimension 1")
#-----------------------------

function foo_par2_codim1(x, y)
    x^6.0 + y^6.0 - 2.0^2.0 +0.0*x*y
    #abs(x)^0.4 + abs((y))^1.4 - 1.0
   # y-abs(sin(x))
    #((x^2.0 + y) - 1.0^2.0)*x #TODO: test with this function, too
end
mymdbm = MDBM_Problem(foo_par2_codim1, [-3.0:3.0, -3.0:3.0])
@time solve!(mymdbm, 2,interpolationorder=1)#number of refinements - increase it slightly to see smoother results 

# mymdbm = MDBM_Problem(foo_par2_codim1, [[-0.0,2.0], [-0.0,2.0]])
# @time solve!(mymdbm, 1,interpolationorder=1)#number of refinements - increase it slightly to see smoother results 




fig = Figure(size = (1300, 1000))  # width x height in pixels
ax1 = GLMakie.Axis(fig[1, 1])
ax2 = GLMakie.Axis(fig[1, 2])
ax3=GLMakie.Axis(fig[1, 3]; limits=((-10.0, 0.0), nothing))

xy_sol = getinterpolatedsolution(mymdbm)
scatter!(ax1,xy_sol...)
nc=mymdbm.ncubes[10]
MDBM.getscaled_local_point(MDBM.ncube_error_vector(nc),nc,mymdbm.axes)


using LinearAlgebra
##
#for _ in 1:5

MDBM.doubling!(mymdbm, [1,2])
nc_list=1:size(mymdbm.ncubes, 1)


errv_s =[MDBM.getscaled_local_point(MDBM.ncube_error_vector(nc),nc,mymdbm.axes)  for nc in mymdbm.ncubes]
err_norm = norm.(errv_s)
nc_list= nc_list[err_norm .> sort(err_norm)[(length(err_norm) *4)÷ 8]]
#nc_list= nc_list[err_norm .> 0.001]
@show length(nc_list)
#MDBM.refinencubes!(mymdbm.ncubes, 1:size(mymdbm.ncubes, 1) ÷ 2, [1, 2])
MDBM.refinencubes!(mymdbm.ncubes,  nc_list, [1,2])
#MDBM.refinencubes!(mymdbm.ncubes, 1:size(mymdbm.ncubes, 1) , [1])
#MDBM.refinencubes!(mymdbm.ncubes, 1:size(mymdbm.ncubes, 1) , [2])
MDBM.interpolate!(mymdbm,interpolationorder=1)
#end


MDBM.getcornerval(mymdbm.ncubes[1], mymdbm)
allcorners = MDBM.corner(mymdbm.ncubes, mymdbm.T01)

x = reduce(vcat, [[[mymdbm.axes[1][xy[1]] for xy in getindex(nc, [1, 2, 4, 3, 1])]..., NaN] for nc in allcorners])
y = reduce(vcat, [[[mymdbm.axes[2][xy[2]] for xy in getindex(nc, [1, 2, 4, 3, 1])]..., NaN] for nc in allcorners])

#f = Figure()
lines!(ax1,x, y)
#lines!(x, y)

# # n-cube interpolation
# xy_sol = getinterpolatedsolution(mymdbm)
# scatter!(xy_sol..., markersize=15, color=:yellow, marker='o', strokewidth=3, label="solution")


MDBM.ncube_error_vector.(mymdbm.ncubes)

errv_s =(MDBM.ncube_error_vector.(mymdbm.ncubes))
err_norm = norm.(errv_s)
xy_sol = getinterpolatedsolution(mymdbm)
ms=((err_norm./maximum(err_norm))).^0.2 .*10

empty!(ax1)
sc = scatter!(ax1,xy_sol..., color =err_norm)#, markersize =ms)
Colorbar(fig[1, 2], sc)


hist!(ax3, log.(err_norm); bins = 30)#, color = :steelblue)

@show fig


#end


##
#--------------------------- Sub-cube interpolation----------------
ax2 = GLMakie.Axis(f[1, 2]; xminorticks=mymdbm.axes[1].ticks, yminorticks=mymdbm.axes[2].ticks, kwargs...)

#calcuatin the sub-cubes interpolations stored in the mymdbm.ncubes[i].posinterp
@time interpsubcubesolution!(mymdbm)
#extracting the resutls to from the 
path2points = extract_paths(mymdbm)

#extracting the unique points and plotting
puniq = unique(collect(Iterators.flatten(Iterators.flatten(path2points))))
scatter!(getindex.(puniq, 1), getindex.(puniq, 2), label="subface - solution")



#exctracing the simplexes for each ncube
flatened_path2points = collect(Iterators.flatten(path2points))
#eliminating the points with less than 2 points (caused by fininte precision)
true_truflatened_path2points = flatened_path2points[length.(flatened_path2points).==2]
#plotting the lines between the points
lines2plot = [(Point2f(ploc[1]), Point2f(ploc[2])) for ploc in true_truflatened_path2points]
linesegments!(lines2plot, label="subface - connection")


display(GLMakie.Screen(), f)

axislegend(ax1)
axislegend(ax2)