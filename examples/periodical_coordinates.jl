
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

#Brute Force plotting for reference
ϕv=LinRange(0.0, 2pi,1000)
ϕv=LinRange(-2pi, 2pi,1000)
θv=LinRange(-pi/2, pi/2,1000)
s=[foo(angle2pos_sphere(ϕ ,θ)...) for ϕ in ϕv,  θ in θv]
xyz=[angle2pos_sphere(ϕ ,θ) for ϕ in ϕv,  θ in θv]


#mymdbm = MDBM_Problem(foo, [-3.0:3.0, -3.0:3.0, -3.0:3.0])
#solve!(mymdbm, 3, verbosity=1) #number of refinements - increase it slightly to see smoother results 

f = Figure(size=(1000, 600))
#ax1 = GLMakie.Axis3(f[1, 1])
contour(f[1,1],ϕv,θv,s,levels=[0.0,0.00000])

#-----------------------------
#MDBM with periodic coordinates

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



#---------------------- Plot the results on the sphere -------------------


xyz_2D_sol=angle2pos_sphere.(xyz_sol...) 
scatter(f[1,2],getindex.(xyz_2D_sol,1),getindex.(xyz_2D_sol,2),getindex.(xyz_2D_sol,3),color=:red)

xyz_2D_eval=angle2pos_sphere.(xy_val...) 
fval=getevaluatedfunctionvalues(mymdbm_SP)
scatter!(f[1,2],getindex.(xyz_2D_eval,1),getindex.(xyz_2D_eval,2),getindex.(xyz_2D_eval,3),color=sign.(fval))
f



