
struct MDBMcontainer{RTf,RTc,AT}
    funval::RTf
    cval::RTc
    callargs::AT
end

#TODO memoization ofr multiple functions Vector{Function}
struct MemF{fT,cT,RTf,RTc,AT} <:Function
    f::fT
    c::cT
    fvalarg::Vector{MDBMcontainer{RTf,RTc,AT}}#funval,callargs
    memoryacc::Vector{Int64} #number of function value call for already evaluated parameters
    MemF(f::fT,c::cT,cont::Vector{MDBMcontainer{RTf,RTc,AT}}) where {fT,cT,RTf,RTc,AT}=new{fT,cT,RTf,RTc,AT}(f,c,cont,[Int64(0)])
end

(memfun::MemF{fT,cT,RTf,RTc,AT})(::Type{RTf},::Type{RTc},args...,) where {fT,cT,RTf,RTc,AT} =( memfun.f(args...,)::RTf, memfun.c(args...,)::RTc)

function (memfun::MemF{fT,cT,RTf,RTc,AT})(args...,) where {fT,cT,RTf,RTc,AT}
    location=searchsortedfirst(memfun.fvalarg,args,lt=(x,y)->isless(x.callargs,y));
    if length(memfun.fvalarg)<location
        x=memfun(RTf,RTc,args...,);
        push!(memfun.fvalarg,MDBMcontainer{RTf,RTc,AT}(x...,args))
        return x
    elseif  memfun.fvalarg[location].callargs!=args
        x=memfun(RTf,RTc,args...,);
        insert!(memfun.fvalarg, location, MDBMcontainer{RTf,RTc,AT}(x...,args))
        return x
    else
        memfun.memoryacc[1]+=1;
        return (memfun.fvalarg[location].funval,memfun.fvalarg[location].cval);
    end
end

function (memfun::MemF{fT,cT,RTf,RTc,AT})(args::Tuple) where {fT,cT,RTf,RTc,AT}
    memfun(args...,)
end

"""
    Axis{T}

Represent the parameter space in a discretized grid (vector) in `ticks`.
"""
struct Axis{T} <: AbstractVector{T}
    ticks::Vector{T}
    name
    function Axis(T::Type,a::AbstractVector,name=:unknown)
        new{T}(T.(a), name)
    end
end
Base.getindex(ax::Axis{T}, ind) where T = ax.ticks[ind]::T
Base.setindex!(ax::Axis, X, inds...) = setindex!(ax.ticks, X, inds...)
Base.size(ax::Axis) = size(ax.ticks)

function Axis(a::AbstractVector{T}, name=:unknown) where T<:Real
    Axis(Float64, a, name)
end
function Axis(a::AbstractVector{T}, name=:unknown) where T
    Axis(T, a, name)
end
function Axis(a, name=:unknown) where T
    Axis([a...], name)
end
function Axis(a::Axis)
    a
end

struct Axes{N,aT}
    axs::aT
end
function Axes(axs...)
    axstosave=Axis.(axs)
    Axes{length(axstosave),typeof(axstosave)}(Axis.(axstosave))
end
(axs::Axes)(ind...) = getindex.(axs.axs,ind)
# (axs::Axes)(ind::AbstractVector{<:Integer}) = axs(ind...)
Base.getindex(axs::Axes,inds...) = axs.axs[inds...]
Base.getindex(axs::Axes,inds::Vector{Integer}) = axs.axs[inds]
Base.iterate(axs::Axes) = iterate(axs.axs)
Base.iterate(axs::Axes,i) = iterate(axs.axs,i)
Base.size(::Axes{N,aT}) where {N,aT} = (N,)
Base.size(::Axes{N,aT},::Integer) where {N,aT} = N
Base.length(::Axes{N,aT}) where {N,aT} = N


struct NCube{IT,FT,N}
    corner::MVector{N,IT} #"bottom-left" #Integer index of the axis
    size::MVector{N,IT}#Integer index of the axis
    posinterp::MVector{N,FT}#relative coordinate within the cube "(-1:1)" range
    bracketingncube::Bool
    # gradient ::MVector{MVector{T}}
    # curvnorm::Vector{T}
end


Base.isless(a::NCube{IT,FT,N},b::NCube{IT,FT,N}) where IT where FT where N = Base.isless([a.corner,a.size],[b.corner,b.size])
Base.isequal(a::NCube{IT,FT,N},b::NCube{IT,FT,N}) where IT where FT where N = all([a.corner==b.corner,a.size==b.size])
import Base.==
==(a::NCube{IT,FT,N},b::NCube{IT,FT,N}) where IT where FT where N = all([a.corner==b.corner,a.size==b.size])


"""
    struct MDBM_Problem{fcT,N,Nf,Nc,t01T,t11T,IT,FT}

Store the main data for the Multi-Dimensional Bisection Method.

# Examples
```julia
include("MDBM.jl")
using Reexport
@reexport using .MDBM

function foo(x,y)
    x^2.0+y^2.0-2.0^2.0
end
function c(x,y) #only the c>0 domain is analysed
    x-y
end

ax1=Axis([-5,-2.5,0,2.5,5],"x") # initial grid in x direction
ax2=Axis(-5:2:5.0,"y") # initial grid in y direction

mymdbm=MDBM_Problem(foo,[ax1,ax2],constraint=c)
```
"""
struct MDBM_Problem{fcT,N,Nf,Nc,t01T,t11T,IT,FT,aT}
    "(memoized) function and constraint in a Tuple"
    fc::fcT
    "The vector of discretized parameters space"
    axes::Axes{N,aT}
    "the bracketing n-cubes (which contains a solution)"
    ncubes::Vector{NCube{IT,FT,N}}
    T01::t01T#
    T11pinv::t11T
end

function MDBM_Problem(fc::fcT,axes,ncubes::Vector{NCube{IT,FT,N}},Nf,Nc) where fcT<:Function where IT<:Integer where FT<:AbstractFloat where N
    T01=T01maker(Val(N))
    T11pinv=T11pinvmaker(Val(N))
    MDBM_Problem{fcT,N,Nf,Nc,typeof(T01),typeof(T11pinv),IT,FT,typeof((axes...,))}(fc,Axes(axes...),
    [NCube{IT,FT,N}(SVector{N,IT}([x...]),SVector{N,IT}(ones(IT,length(x))),SVector{N,FT}(zeros(IT,length(x))),true) for x in Iterators.product((x->1:(length(x.ticks)-1)).(axes)...,)][:]
    ,T01,T11pinv)
end

function MDBM_Problem(f::Function, axes0::AbstractVector{<:Axis};constraint::Function=(x...,)->nothing, memoization::Bool=true,
    #Nf=length(f(getindex.(axes0,1)...)),
    Nf=f(getindex.(axes0,1)...) === nothing ? 0 : length(f(getindex.(axes0,1)...)),
    Nc=constraint(getindex.(axes0,1)...) === nothing ? 0 : length(constraint(getindex.(axes0,1)...)))#Float16(1.), nothing
    axes=deepcopy.(axes0);
    argtypesofmyfunc=map(x->typeof(x).parameters[1], axes);#Argument Type
    AT=Tuple{argtypesofmyfunc...};
    type_f=Base.return_types(f,AT)
    if length(type_f)==0
        error("input of the function is not compatible with the provided axes")
    else
        RTf=type_f[1];#Return Type of f
    end

    type_con=Base.return_types(constraint, AT)
    if length(type_con)==0
        error("input of the constraint function is not compatible with the provided axes")
    else
        RTc=type_con[1];#Return Type of the constraint function
    end

    if memoization
        fun=MemF(f, constraint,Array{MDBMcontainer{RTf,RTc,AT}}(undef, 0));
    else
        fun=(x)->(f(x...), constraint(x...));
    end
    Ndim=length(axes)
    MDBM_Problem(fun, axes, Vector{NCube{Int64,Float64,Ndim}}(undef, 0),  Nf, Nc)
end

function MDBM_Problem(f::Function, a::AbstractVector{<:AbstractVector};constraint::Function=(x...,)->nothing, memoization::Bool=true,
    #Nf=length(f(getindex.(axes0,1)...)),
    Nf=f(getindex.(a,1)...) === nothing ? 0 : length(f(getindex.(a,1)...)),
    Nc=constraint(getindex.(a,1)...) === nothing ? 0 : length(constraint(getindex.(a,1)...)))
    axes=[Axis(ax) for ax in a]
    MDBM_Problem(f, axes, constraint = onstraint, memoization=memoization, Nf=Nf, Nc=Nc)#,Vector{NCube{Int64,Float64,Val(Ndim)}}(undef, 0))
end


function Base.show(io::IO, mdbm::MDBM_Problem{fcT,N,Nf,Nc,t01T,t11T,IT,FT,aT}) where {fcT,N,Nf,Nc,t01T,t11T,IT,FT,aT}
    println(io,"Multi-Dimensional Bisection Method Problem")
    println(io,"  parameter dimension:   ", N)#typeof(mdbm).parameters[2])#N
    println(io,"  co-dimension:          ", Nf)#typeof(mdbm).parameters[3])#Nf
    println(io,"  number of constraints: ", Nc)#typeof(mdbm).parameters[4])#Nc
    println(io,"Axes:")
    for k in 1:length(mdbm.axes)
        println(io,"  axis #",k,": elements: ",length(mdbm.axes[k]),"; elementtype: ", typeof(mdbm.axes[k][1]),"; first: ", mdbm.axes[k][1],"; last: ",mdbm.axes[k][end])
    end
        println(io,"number of bracketing n-cubes: ", length(mdbm.ncubes))
        println()

    if (typeof(mdbm.fc) <: MemF)
        println(io,"number of function evaluation: ", length(mdbm.fc.fvalarg))
        println(io,"number of memoized function call: ", mdbm.fc.memoryacc[1])
        #println(io,"Ration of function evaluation compared to 'brute-force method': ", length(mdbm.fc.fvalarg)/prod(length.(mdbm.axes)))
    else
     println(io,"non-memoized version")
    end
end
