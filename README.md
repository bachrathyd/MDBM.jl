# MDBM.jl

[![Powered by NumFOCUS](https://img.shields.io/badge/powered%20by-NumFOCUS-orange.svg?style=flat&colorA=E1523D&colorB=007D8A)](http://numfocus.org)

Multi-Dimensional Bisection Method (MDBM) is an efficient and robust root-finding algorithm, which can be used to determine whole high-dimensional submanifolds (points, curves, surfaces…) of the roots of implicit non-linear equation systems, even in cases, where the number of unknowns surpasses the number of equations.

<a href="https://www.codecogs.com/eqnedit.php?latex=f_i(x_j)=0&space;\quad&space;i=1...k&space;\quad&space;j=1...l,&space;\quad&space;k&space;\leq&space;l" target="_blank"><img src="https://latex.codecogs.com/gif.latex?f_i(x_j)=0&space;\quad&space;i=1...k&space;\quad&space;j=1...l,&space;\quad&space;k&space;\leq&space;l" title="f_i(x_j)=0 \quad i=1...k \quad j=1...l, \quad k \leq l" /></a>

This type of problems can be found in many different field of science, just to mention a few:
- differential geometry (isolines, isosurfaces in higher dimensions)
- analysis of linkages (mechanical: workspace of robots)
- stability computation, stabilyzability diagrams


This method is an alternative to the contour plot or to isosurfaces in higher dimension, however, it has as the advantage of being able to handle multiple functions at once. <br>
In addition, it uses far less function evaluation than the brute-force approach, making it much faster and more memory efficient, especially for complex tasks.


## Introduction

The bisection method - or the so-called interval halving method - is one of the simplest root-finding algorithms which is used to find zeros of continuous non-linear functions.
This method is very robust and it always tends to the solution if the signs of the function values are different at the borders of the chosen initial interval.

Geometrically, root-finding algorithms of __f__(__x__)=0 find one intersection point of the graph of the function with the axis of the independent variable.
In many applications, this 1-dimensional intersection problem must be extended to higher dimensions, e.g.: intersections of surfaces in a 3D space (volume), which can be described as a system on non-linear implicit equations. In higher dimensions, the existence of multiple solutions becomes very important, since the intersection of two surfaces can create multiple intersection curves.

MDBM algorithm canhandle automatically:
- multiple solutions 
- arbitrary number of parameter (typically: 3-6)
- arbitrary number implicit equations
- multiple constraints in the parameter space
- handle degenerated functions
- first order interpolation (and convergence rate)
- provides the gradients of the equations for the roots


## Citing
The software in this ecosystem was developed as part of academic research. If you use the MDBM.jl package as part of your research, teaching, or other work, I would be grateful if you could cite my corresponding publication: <https://pp.bme.hu/me/article/view/1236/640>


## Web:
<https://www.mm.bme.hu/~bachrathy/research_EN.html>

# Quick start

Preparation
```julia
include("MDBM.jl")
using Reexport
@reexport using .MDBM

using PyPlot;
pygui(true);
```
## Example 1
Computation of a circle section defined by the implicit equation `foo(x,y) == 0` and by the constraint `c(x,y) > 0`
```julia
function foo(x,y)
    x^2.0+y^2.0-2.0^2.0
end
function c(x,y) #only the c>0 domain is analysed
    x-y
end

ax1=Axis([-5,-2.5,0,2.5,5],"x") # initial grid in x direction
ax2=Axis(-5:2:5.0,"y") # initial grid in y direction

mymdbm=MDBM_Problem(foo,[ax1,ax2],constraint=c)
iteration=5 #number of refinements (resolution doubling)
solve!(mymdbm,iteration) 


#points where the function foo was evaluated
x_eval,y_eval=getevaluatedpoints(mymdbm)

#interpolated points of the solution (approximately where foo(x,y) == 0 and c(x,y)>0)
x_sol,y_sol=getinterpolatedsolution(mymdbm)

fig = figure(1);clf()
scatter(x_eval,y_eval,s=5)
scatter(x_sol,y_sol,s=5)
```
<img src="assets/circe_2D_points.png"
     alt="solution "/>

Perform the line connection
```julia
myDT1=connect(mymdbm);
for i in 1:length(myDT1)
    dt=myDT1[i]
    P1=getinterpolatedsolution(mymdbm.ncubes[dt[1]],mymdbm)
    P2=getinterpolatedsolution(mymdbm.ncubes[dt[2]],mymdbm)
    plot([P1[1],P2[1]],[P1[2],P2[2]], color="k")
end 
```

<img src="assets/circe_2D_line.png"
     alt="solution "/>

## Example 2
Intersection of two sphere in 3D
```julia
using LinearAlgebra

axes=[-2:2,-2:2,-2:2]

fig = figure(3);clf()


#Sphere1
fS1(x...) = norm(x,2.0)-1

Sphere1mdbm=MDBM_Problem(fS1,axes)
solve!(Sphere1mdbm,4)
a_sol,b_sol,c_sol=getinterpolatedsolution(Sphere1mdbm)
plot3D(a_sol,b_sol,c_sol,linestyle="", marker=".", markersize=1);


#Sphere2
fS2(x...) = norm(x .- 0.5, 2.0) -1.0

Sphere2mdbm=MDBM_Problem(fS2,axes)
solve!(Sphere2mdbm,4)
a_sol,b_sol,c_sol=getinterpolatedsolution(Sphere2mdbm)
plot3D(a_sol,b_sol,c_sol,linestyle="", marker=".", markersize=1);

#Intersection
fS12(x...) = (fS1(x...), fS2(x...))

Intersectmdbm=MDBM_Problem(fS12,axes)
solve!(Intersectmdbm,6)
a_sol,b_sol,c_sol=getinterpolatedsolution(Intersectmdbm)
plot3D(a_sol,b_sol,c_sol,color="k",linestyle="", marker=".", markersize=2);

```

<img src="assets/sphere_intersection.png"
     alt="solution "/>

## Example 3
Mandelbrot set (example for a non-smooth problem)

```julia
function mandelbrot(x,y)    
    c=x+y*im
    z=Complex(zero(c))
    k=0
    maxiteration=1000
    while (k<maxiteration && abs(z)<4)
            z=z^2+c
            k=k+1
        end
    return abs(z)-2
end

Mandelbrotmdbm=MDBM_Problem(mandelbrot,[-5:2,-2:2])
solve!(Mandelbrotmdbm,9)
a_sol,b_sol=getinterpolatedsolution(Mandelbrotmdbm)
fig = figure(3);clf()
plot(a_sol,b_sol,linestyle="", marker=".", markersize=1)
```
<img src="assets/Mandelbrot.png"
     alt="solution "/>

## History

I am an assistant professor at the Budapest University of Technology and Economics, at the Faculty of Mechanical Engineering, the Department of Applied Mechanics.
During my studies and research, I have to determine stability charts of models described by delayed differential equations, which are typically formed as a "3 parameter / 2 implicit equation" problem. I have faced the difficulty that there is no applicable solution in any available software (e.g.: Mathematica, Matlab,...) which could easily be used in engineering problems. 
Due to this reason, I have started to develop the Multi-Dimensional Bisection Method since 2006 in Matlab, and I have been improving it since, by adding new features from time to time.  
I hope that others also find this package a useful tool in their work.

Best regards,
Dr. Daniel Bachrathy

#

![NumFOCUS logo](assets/numfocus-logo.png)

MDBM is a fiscally sponsored project of [NumFOCUS](https://numfocus.org), a
nonprofit dedicated to supporting the open source scientific computing
community.