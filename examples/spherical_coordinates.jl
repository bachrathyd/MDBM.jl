
using MDBM
using GLMakie

GLMakie.closeall()
GLMakie.activate!(;title = "2 parameters, codimension 1")
#-----------------------------

function angle2pos_sphere(ϕ ,θ)
    dfi=0.5
    x=cos(ϕ+dfi)*cos(θ)
    y=sin(ϕ+dfi)*cos(θ)
    z=sin(θ)
    return x,y,z
end
function foo(x,y,z)
    p=0.5
    abs(x-0.75)^p+abs(y-0.7)^p+abs(z-0.4)^p-0.85^p
end

ϕv=LinRange(0.0, 2pi,1000)
ϕv=LinRange(-2pi, 2pi,1000)
θv=LinRange(-pi/2, pi/2,1000)
s=[foo(angle2pos_sphere(ϕ ,θ)...) for ϕ in ϕv,  θ in θv]
xyz=[angle2pos_sphere(ϕ ,θ) for ϕ in ϕv,  θ in θv]


#mymdbm = MDBM_Problem(foo, [-3.0:3.0, -3.0:3.0, -3.0:3.0])
#solve!(mymdbm, 3, verbosity=1) #number of refinements - increase it slightly to see smoother results 

f = Figure(size=(1000, 600))
#ax1 = GLMakie.Axis3(f[1, 1])
contour(f[1,1],ϕv,θv,s,levels=[0.0,0.1])



ϕv=MDBM.Axis(LinRange(0.0, 2pi,20),"ϕ",true)
#ϕv=MDBM.Axis(LinRange(2pi, 0.0,20),"ϕ",false)

θv=MDBM.Axis(LinRange(-pi/2, pi/2,8))
Foo_sphere(x...)=foo(angle2pos_sphere(x...)...)
Foo_sphere(1.2,3.2)
Foo_sphere(0.4,0.2)
mymdbm_SP = MDBM_Problem(  Foo_sphere, [ϕv,θv])
solve!(mymdbm_SP, 8, verbosity=1) #number of refinements - increase it slightly to see smoother results 
# n-cube interpolation
xyz_sol = getinterpolatedsolution(mymdbm_SP)
scatter!( xyz_sol..., markersize=6, color=:red, marker='x', strokewidth=3, label="solution")


xy_val = getevaluatedpoints(mymdbm_SP)
fval=getevaluatedfunctionvalues(mymdbm_SP)
scatter!(xy_val...,color=sign.(fval),label = "evaluated")

#----------------------

xyz_2D_sol=angle2pos_sphere.(xyz_sol...) 
scatter(f[1,2],getindex.(xyz_2D_sol,1),getindex.(xyz_2D_sol,2),getindex.(xyz_2D_sol,3),color=:red)

xyz_2D_eval=angle2pos_sphere.(xy_val...) 
fval=getevaluatedfunctionvalues(mymdbm_SP)
scatter!(f[1,2],getindex.(xyz_2D_eval,1),getindex.(xyz_2D_eval,2),getindex.(xyz_2D_eval,3),color=sign.(fval))
f








































mymdbm = MDBM_Problem(foo_par2_codim1, [-3.1:3.0, -3.1:3.0])
@time solve!(mymdbm, 5)#number of refinements - increase it slightly to see smoother results 

f = Figure()
#show the final resolution of the grid based on the minorticks
kwargs = (; xminorticksvisible = true, xminorgridvisible = true, yminorticksvisible = true, yminorgridvisible = true)
ax1=GLMakie.Axis(f[1, 1]; xminorticks = mymdbm.axes[1].ticks, yminorticks  = mymdbm.axes[2].ticks, kwargs...)

# n-cube interpolation
xy_sol = getinterpolatedsolution(mymdbm)
scatter!(xy_sol..., markersize = 15, color = :red,marker ='x',strokewidth=3,label = "solution")

# show the points where the function is evaluated
xy_val = getevaluatedpoints(mymdbm)
fval=getevaluatedfunctionvalues(mymdbm)
scatter!(xy_val...,color=sign.(fval),label = "evaluated")

# connecting and plotting the "mindpoints" of the n-cubes
DT1 = connect(mymdbm)
edge2plot_xyz = [reduce(hcat, [i_sol[getindex.(DT1, 1)], i_sol[getindex.(DT1, 2)], fill(NaN, length(DT1))])'[:] for i_sol in xy_sol]
lines!(edge2plot_xyz..., linewidth=5,label = "midpoints solution connected")


#plotting the gradintes
gxyz=getinterpolatedgradient(mymdbm.ncubes,mymdbm)
arrows!(xy_sol..., gxyz[1]..., arrowsize = 0.01, lengthscale = 0.1,label = "gradient")#    arrowcolor = strength, linecolor = strength)



#--------------------------- Sub-cube interpolation----------------
ax2=GLMakie.Axis(f[1, 2]; xminorticks = mymdbm.axes[1].ticks, yminorticks  = mymdbm.axes[2].ticks, kwargs...)

#calcuatin the sub-cubes interpolations stored in the mymdbm.ncubes[i].posinterp
interpsubcubesolution!(mymdbm)
#extracting the resutls to from the 
path2points = extract_paths(mymdbm)

#extracting the unique points and plotting
puniq = unique(collect(Iterators.flatten(Iterators.flatten(path2points))))
scatter!(getindex.(puniq, 1), getindex.(puniq, 2),label = "subface - solution")



#exctracing the simplexes for each ncube
flatened_path2points = collect(Iterators.flatten(path2points))
#eliminating the points with less than 2 points (caused by fininte precision)
true_truflatened_path2points = flatened_path2points[length.(flatened_path2points) .== 2]
#plotting the lines between the points
lines2plot = [(Point2f(ploc[1]) , Point2f(ploc[2])) for ploc in true_truflatened_path2points]
linesegments!(lines2plot,label = "subface - connection")


display(GLMakie.Screen(), f)

axislegend(ax1)
axislegend(ax2)