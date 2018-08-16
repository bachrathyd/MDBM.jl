struct Axis <: AbstractVector{AbstractFloat}
    ticks::Vector{<:AbstractFloat}
    name::Symbol
end
Base.getindex(ax::Axis, ind) = ax.ticks[ind]
Base.setindex!(ax::Axis, X, inds...) = setindex!(ax.ticks, X, inds...)
Base.size(ax::Axis) = size(ax.ticks)
function Axis(a::Vector{<:AbstractFloat})
    Axis(a, :unknown)
end
function Axis(a::AbstractVector, name::Symbol=:unknown;dtype::Type=Float64)
    Axis(convert(Vector{dtype}, a), name)
end
function Axis(a::Axis)
    a
end
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#

struct Constraint
    g::Function
end
Constraint(c::Constraint)=c
(c::Constraint)(t) = c.g(t)
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#

struct Problem
    f::Function
    c::Constraint
    axes::Vector{Axis}
    is_constrained::Bool
end
function Problem(f::Function, axes::Vector{<:AbstractVector})
    Problem(f, Constraint(x -> error("NO CONSTRAINT DEFINED")), Axis.(axes), false)
end
function Problem(f::Function,g::Constraint, axes::Vector{<:AbstractVector})
    Problem(f, Constraint(g), Axis.(axes), true)
end

#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#

# mutable struct Results (??????)
struct Results
    problem::Problem
    #...
end

mutable struct mdbm_object
    f::Function## evaluated function
    fconstrain::Function## evaluated function
    fvectorized::Bool ## is function f can be called in a vectorized form
    ax::Array{Array{Float64,1},1}##description of the initial grid
    Ndim::Int64
    Nax::Array{Int64,1}
    Naxstride::Array{Int64,1}
  
    Ncodim::Int64
    HC::Array{Float64,2}#T computed function values and constrain value
    linindex::Array{Int64,1} # corresponding linear-indices
    vectindex::Array{Int64,2}# corresponding sub-indices
    #N::Int64
    #compind
    #pointerp
    #DT
    ncubelin::Array{Int64,1} # linear-indices of the bracketing n-cubes
    ncubevect::Array{Int64,2} # sub-indices of the bracketing n-cubes
  
    posinterp::Array{Float64,2}# the interpolated valuse within the bracketing n-cubes
    gradient::Array{Float64,3}# the corresponding gradients within the bracketing n-cubes
  
    DT1::Array{Int64,2}#line 'tiangulation' of the resultant interpolated values (basd on the n-cube sub-indices)
    DT2::Array{Int64,2}#surface tiangulation of the resultant interpolated values (basd on the n-cube sub-indices)
  
    isconstrained::Bool #is f provides contrain? all the constraints are combinded!!! length C===1
    interporder::Int64 # 0,1,2
    selectionmode::Int64 # 0-safe selection, 1st order interpolation mased,2???
end