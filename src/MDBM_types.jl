
struct MemF0{FT,CT}
    f::Function
    funval::Vector{FT}
    callargs::Vector{CT}
    #MemF(fun::Function)= MemF(fun::Function,Any[],Any[])
end


function (memfun::MemF0)(args...)
    println(args...)
    println(memfun.callargs)
    #location=findfirst(x -> x==args...,memfun.callargs);
    location=findfirst(x -> x==args,memfun.callargs);
    println(location)
    if location==nothing
        println("ujraszamolas")
        x=memfun.f(args...);
        println(x)
        push!(memfun.funval,x)
        push!(memfun.callargs,args);
        return x
    else
        return memfun.funval[location]
    end
end



struct MemF{FT}  #<:Function
    f::Function
    fvalarg::Vector{FT}#funval,callargs
    #MemF(fun::Function)= MemF(fun::Function,Any[],Any[])
end



NN=200_000
N=200


bb=x->sin(x)
Base.return_types(bb,(Int32,))

function (memfun::MemF)(args...)#non sorted searching
    location=findfirst(x -> x[2]==args,memfun.fvalarg);
    if location==nothing
        x=memfun.f(args...);
        push!(memfun.fvalarg,(x,args))
        #sort!(memfun.fvalarg,by=x -> x[2]);
        return x
    else
        return memfun.fvalarg[location][1]
    end
end



fa2=MemF((x,y)->(x^2.0)/y,Array{Tuple{Float64,Tuple{Float16, Float16}}}(undef, 0))

fa2(-1,-1.)
@time for k in 1:NN
    fa2(ceil(rand(Float16)*N),ceil(rand(Float16)*N))
end
println(length(fa2.fvalarg))


function (memfun::MemF)(args...)
    #location=findfirst(x -> x[2]==args,memfun.fvalarg);
    #location=searchsortedfirst([x[2] for x in memfun.fvalarg],args);
    location=searchsortedfirst(memfun.fvalarg,args,lt=(x,y)->isless(x[2],y));
    #println(location)
    #if location==nothing
    if length(memfun.fvalarg)<location
        #println("ujraszamolas a végére")
        x=memfun.f(args...);
        push!(memfun.fvalarg,(x,args))
        return x
    elseif  memfun.fvalarg[location][2]!=args
        #println("ujraszamolas közé illeszt")
        x=memfun.f(args...);

        # splice!(memfun.fvalarg,location,[(x,args),memfun.fvalarg[location]])
        # TODO Insert
        # insert!();%%%%%%%%%%%%
        insert!(memfun.fvalarg, location, (x,args))
        return x
    else
        #println("mar megvolt")
        return memfun.fvalarg[location][1]
    end
end
#        insert_and_dedup!(v::Vector, x) = (splice!(v, searchsorted(v,x), [x]); v)

fa3=MemF((x,y)->(x^2.0)/y,Array{Tuple{Float64,Tuple{Float16, Float16}}}(undef, 0))

fa3(-1,-1.)
@time for k in 1:NN
    fa3(ceil(rand(Float16)*N),ceil(rand(Float16)*N))
end
println(length(fa3.fvalarg))


@time for k in 1:NN
    fa3.f(ceil(rand(Float16)*N),ceil(rand(Float16)*N))
end
println("----------------")


function bbb(x,y)
    (x^2)*y
end
bbb(4.3,2)
Base.return_types(fa3,(Int32,Float64))






struct Myfun
    f::Function
end

function Base.getindex(mfun::Myfun, ind)
     mfun.f(ind...)
end


function (x::Myfun)(y::Int64)
    println(y)
end



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

struct EVP{restT,axT}
    ax::Vector{axT}
    H::Vector{restT}
    C::Vector{restT}
end

struct mdbm_object{T}
    f::Function## evaluated function
    fconstrain::Function## evaluated function
    fvectorized::Bool ## is function f can be called in a vectorized form
    ax::Array{Array{Float64,1},1}##description of the initial grid
    Ndim::Int64
    Nax::Array{Int64,1}
    Naxstride::Array{Int64,1}

    myevp::Vector{EVP}

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
