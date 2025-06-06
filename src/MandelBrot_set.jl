using MDBM
using GLMakie

GLMakie.closeall()
GLMakie.activate!(;title = "2 parameters, codimension 1 - Mandelbrot set")
#-----------------------------

function mandelbrot(x,y)    
    c = x + y*im
    z = zero(c)
    k = 0
    maxiteration = 1000
    while (k < maxiteration && abs(z) < 4)
        z = z^2 + c
        k += 1
    end
    return abs(z) - 2
end
Mandelbrotmdbm = MDBM_Problem(mandelbrot, [-2:0.5:1, -1.5:0.5:1.5])
solve!(Mandelbrotmdbm, 8)

checkneighbour!(Mandelbrotmdbm,interpolationorder=1,normp=Inf, ncubetolerance=1.9)
interpolate!(Mandelbrotmdbm,interpolationorder=1,normp=Inf, ncubetolerance=0.8)

f = Figure(size=(1500,800))
#show the final resolution of the grid based on the minorticks
kwargs = (; xminorticksvisible = true, xminorgridvisible = true, yminorticksvisible = true, yminorgridvisible = true)
ax1=GLMakie.Axis(f[1, 1]; xminorticks = Mandelbrotmdbm.axes[1].ticks, yminorticks  = Mandelbrotmdbm.axes[2].ticks, kwargs...)

# n-cube interpolation
xy_sol = getinterpolatedsolution(Mandelbrotmdbm)

scatter!(xy_sol..., markersize = 10, color = :red,marker ='x',label = "solution")


# show the points where the function is evaluated
xy_val = getevaluatedpoints(Mandelbrotmdbm)
fval=getevaluatedfunctionvalues(Mandelbrotmdbm)
scatter!(xy_val...,color=sign.(fval),label = "evaluated")


display(GLMakie.Screen(), f)

