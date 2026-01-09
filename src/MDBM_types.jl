"""
    SortedCache{K,V}

A small wrapper that keeps:
  • `data::Dict{K,V}` for O(1) lookup/insert  
  • `keys::Vector{K}` (sorted) for position‐indexing/searchsortedfirst  

Insertion only does a binary‐search+`insert!` on the *key* array (cheap `K` copies,
not shifting huge containers), and storing the heavy `V` in a hash‐table.
"""
struct SortedCache{K,V}
    data::Dict{K,V}
    keys::Vector{K}
end

SortedCache{K,V}() where {K,V} = SortedCache(Dict{K,V}(), Vector{K}())

# insert or update
function insert!(sc::SortedCache{K,V}, key::K, val::V) where {K,V}
    if !haskey(sc.data, key)
        i = searchsortedfirst(sc.keys, key)
        Base.insert!(sc.keys, i, key)      # shifts only K’s, not whole containers
    end
    sc.data[key] = val
    return nothing
end


"""
    append_SortedCache!(sc::SortedCache{K,V}, kvs::Vector{Pair{K,V}})

Insert many (key ⇒ val) pairs at once.  New keys are sorted,
duplicates overwrite, and the merge cost is O(n₁+n₂), not O(n₂·n₁).
"""
function append_SortedCache!(sc::SortedCache{K,V}, kvs::Vector{Pair{K,V}}) where {K,V}
    # 1) sort your incoming batch by key
    sort!(kvs, by=p -> p.first)

    # 2) extract sorted keys & vals
    newkeys = [p.first for p in kvs]
    newvals = [p.second for p in kvs]

    # 3) merge old sc.keys and newkeys into one vector
    oldkeys = sc.keys
    merged = Vector{K}(undef, length(oldkeys) + length(newkeys))
    i = j = k = 1
    while i ≤ length(oldkeys) && j ≤ length(newkeys)
        if oldkeys[i] === newkeys[j] || oldkeys[i] < newkeys[j]
            merged[k] = oldkeys[i]
            i += 1
        else
            merged[k] = newkeys[j]
            sc.data[newkeys[j]] = newvals[j]  # insert into Dict
            j += 1
        end
        k += 1
    end
    # finish leftovers
    while i ≤ length(oldkeys)
        merged[k] = oldkeys[i]
        i += 1
        k += 1
    end
    while j ≤ length(newkeys)
        merged[k] = newkeys[j]
        sc.data[newkeys[j]] = newvals[j]
        j += 1
        k += 1
    end

    resize!(merged, k - 1)
    # instead of sc.keys = merged, do:
    empty!(sc.keys)               # drop all old entries
    Base.append!(sc.keys, merged)      # fill in the merged key‐list
    return nothing
end


# lookup by key
get(sc::SortedCache{K,V}, key::K, default) where {K,V} = get(sc.data, key, default)
Base.getindex(sc::SortedCache{K,V}, key::K) where {K,V} = sc.data[key]

# “vector” interface: idx→ (funval,cval,callargs)
Base.length(sc::SortedCache) = length(sc.keys)
function Base.getindex(sc::SortedCache{K,V}, i::Integer) where {K,V}
    k = sc.keys[i]
    return sc.data[k]
end

# allow `searchsortedfirst(sc, key)` and iteration
Base.iterate(sc::SortedCache) = iterate(sc.keys)
Base.iterate(sc::SortedCache, st) = iterate(sc.keys, st)
Base.searchsortedfirst(sc::SortedCache{K,V}, key::K) where {K,V} =
    searchsortedfirst(sc.keys, key)











struct MDBMcontainer{RTf,RTc,AT}
    funval::RTf
    cval::RTc
    callargs::AT
end

Base.isless(a::AT, b::MDBMcontainer{RTf,RTc,AT}) where {RTf,RTc,AT} = Base.isless(a, b.callargs)
Base.isless(a::MDBMcontainer{RTf,RTc,AT}, b::AT) where {RTf,RTc,AT} = Base.isless(a.callargs, b)
Base.isless(a::MDBMcontainer{RTf,RTc,AT}, b::MDBMcontainer{RTf,RTc,AT}) where {RTf,RTc,AT} = Base.isless(a.callargs, b.callargs)


#TODO: Shifting leads to memory movements! Fix: switch to an OrderedDict{AT,Tuple{RTf,RTc}} (or even a plain Dict plus a separate sorted key‐list only when you need ordering). Lookups/inserts become amortized O(1), no shifting. 
#TODO memoization for multiple functions Vector{Function}
using FunctionWrappers: FunctionWrapper

struct MemF{RTf,RTc,AT} <: Function
    f::FunctionWrapper{RTf,AT}
    c::FunctionWrapper{RTc,AT}
    fvalarg::SortedCache{AT,Tuple{RTf,RTc}}#Vector{MDBMcontainer{RTf,RTc,AT}}#funval,callargs
    memoryacc::Ref{Int64}#MVector{1,Int64} #number of function value call for already evaluated parameters
    MemF(f, c, cont::SortedCache{AT,Tuple{RTf,RTc}}, ::Type{RTf}, ::Type{RTc}, ::Type{AT}) where {RTf,RTc,AT} =
        new{RTf,RTc,AT}(FunctionWrapper{RTf,AT}(f), FunctionWrapper{RTc,AT}(c), cont, Ref(Int64(0)))
end

(memfun::MemF{RTf,RTc,AT})(::Type{RTf}, ::Type{RTc}, args::AT) where {RTf,RTc,AT} = (memfun.f(args...,)::RTf, memfun.c(args...,)::RTc)::Tuple{RTf,RTc}

function (memfun::MemF{RTf,RTc,AT})(args::AT) where {RTf,RTc,AT}
    if haskey(memfun.fvalarg.data, args)
        memfun.memoryacc[] += 1
        return memfun.fvalarg[args]
    else
        x = memfun(RTf, RTc, args)
        insert!(memfun.fvalarg, args, x)
        return x
    end
    #  # deprecated #  #  @show location = searchsortedfirst(memfun.fvalarg, args)
    #  # deprecated #  #  @show length(memfun.fvalarg)
    #  # deprecated #  #  if length(memfun.fvalarg) < location
    #  # deprecated #  #      #@show args
    #  # deprecated #  #      #@show memfun.fvalarg
    #  # deprecated #  #      #error("this should not happen do to the precalcuation!")
    #  # deprecated #  #      x = memfun(RTf, RTc, args)
    #  # deprecated #  #      #push!(memfun.fvalarg, MDBMcontainer{RTf,RTc,AT}(x..., args))
    #  # deprecated #  #      insert!(memfun.fvalarg, args, x)
    #  # deprecated #  #      return x
    #  # deprecated #  #      #elseif memfun.fvalarg[location].callargs != args
    #  # deprecated #  #      #    #@show args
    #  # deprecated #  #      #    #@show memfun.fvalarg[location]
    #  # deprecated #  #      #    #error("this should not happen do to the precalcuation!")
    #  # deprecated #  #      #    x = memfun(RTf, RTc, args)
    #  # deprecated #  #      #    insert!(memfun.fvalarg, args, x)
    #  # deprecated #  #      #    #insert!(memfun.fvalarg, location, MDBMcontainer{RTf,RTc,AT}(x..., args))
    #  # deprecated #  #      #    return x
    #  # deprecated #  #  else
    #  # deprecated #  #      memfun.memoryacc[] += 1
    #  # deprecated #  #      return (memfun.fvalarg[location].funval, memfun.fvalarg[location].cval)
    #  # deprecated #  #  end
end

function (memfun::MemF{RTf,RTc,AT})(args...) where {RTf,RTc,AT}
    memfun(AT(args))
end

#TODO Inconsistent return in the “batch” call - it should return a vector of the results of the function should be renamed e.g.: prefill!
function (memfun::MemF{RTf,RTc,AT})(Vargs::AbstractVector{AT}) where {RTf,RTc,AT}
    #println("Threaded")
    Vargs = unique(sort(Vargs))
    VargsIndex2compute = .!is_sorted_in_sorted(Vargs, memfun.fvalarg)

    #@show sum(VargsIndex2compute)
    TheContainer_sorted = Array{Tuple{RTf,RTc}}(undef, sum(VargsIndex2compute))
    #Threads.@threads    for (index,args) in enumerate(Vargs[VargsIndex2compute])


    Vargs2compute_sorted = sort(Vargs[VargsIndex2compute])
    Threads.@threads for index in eachindex(Vargs2compute_sorted)
        # @show Vargs2compute[index]
        TheContainer_sorted[index] = memfun(RTf, RTc, Vargs2compute_sorted[index])
    end

    # sort them in-place by key:
    append_SortedCache!(memfun.fvalarg, Pair.(Vargs2compute_sorted, TheContainer_sorted))

    # merge_sorted(memfun.fvalarg, TheContainer) # tooo, slow. Maybe due to the memory allocation and copying
    return nothing
end




# depraceted
#function (memfun::MemF{fT,cT,RTf,RTc,AT})(args::Tuple) where {fT,cT,RTf,RTc,AT}
#    println("----sdfsdfg-----------") 
#       memfun(args...,)
#end

"""
    Axis{T}

Represent the parameter space in a discretized grid (vector) in `ticks`.
"""
struct Axis{T} <: AbstractVector{T}
    ticks::Vector{T}
    name
    periodic::Bool
    function Axis(T::Type, a::AbstractVector, name=:unknown, periodic=false)
        new{T}(T.(a), name, periodic)
    end
end
Base.getindex(ax::Axis{T}, ind) where {T} = ax.ticks[ind]::T
Base.setindex!(ax::Axis, X, inds...) = setindex!(ax.ticks, X, inds...)
Base.size(ax::Axis) = size(ax.ticks)

function Axis(a::AbstractVector{T}, name=:unknown, periodic=false) where {T<:Real}
    Axis(T, a, name, periodic)
end
function Axis(a::AbstractVector{T}, name=:unknown, periodic=false) where {T}
    Axis(T, a, name)
end
function Axis(a, name=:unknown)
    Axis([a...], name)
end
function Axis(a::Axis)
    a
end

struct Axes{N,aT}
    axs::aT
end
function Axes(axs...)
    axstosave = Axis.(axs)
    Axes{length(axstosave),typeof(axstosave)}(Axis.(axstosave))
end
(axs::Axes)(ind...) = getindex.(axs.axs, ind)
# (axs::Axes)(ind::AbstractVector{<:Integer}) = axs(ind...)
Base.getindex(axs::Axes, inds...) = axs.axs[inds...]
Base.getindex(axs::Axes, inds::Vector{Integer}) = axs.axs[inds]
Base.iterate(axs::Axes) = iterate(axs.axs)
Base.iterate(axs::Axes, i) = iterate(axs.axs, i)
Base.size(::Axes{N,aT}) where {N,aT} = (N,)
Base.size(::Axes{N,aT}, ::Integer) where {N,aT} = N
Base.length(::Axes{N,aT}) where {N,aT} = N



struct PositionTree{N,FT}
    p::MVector{N,FT}
    subpoints::Vector{PositionTree{N,FT}}
end
# Constructor for a single position (empty subpoints)
PositionTree(p::MVector{N,T}) where {N,T} = PositionTree{N,T}(p, PositionTree{N,T}[])
# Constructor accepting any AbstractArray, automatically converting to MVector
PositionTree(p::AbstractArray{T}) where {T} = PositionTree(MVector{length(p),T}(p))




struct NCube{IT,FT,N,Nfc}
    corner::MVector{N,IT} #"bottom-left" #Integer index of the axis
    size::MVector{N,IT}#Integer index of the axis
    posinterp::PositionTree{N,FT}#relative coordinate within the cube "(-1:1)" range
    bracketingncube::Bool
    gradient::MMatrix{N,Nfc,FT}#relative coordinate within the cube "(-1:1)" range
    #gradient ::MVector{Nfc,MVector{N,FT}}
    # curvnorm::Vector{T}
end


# Base.isless(a::NCube{IT,FT,N,Nfc}, b::NCube{IT,FT,N,Nfc}) where {IT,FT,N,Nfc} = Base.isless([a.corner, a.size], [b.corner, b.size])
# Base.isequal(a::NCube{IT,FT,N,Nfc}, b::NCube{IT,FT,N,Nfc}) where {IT,FT,N,Nfc} = all([a.corner == b.corner, a.size == b.size])
# import Base.==
# ==(a::NCube{IT,FT,N,Nfc}, b::NCube{IT,FT,N,Nfc}) where {IT,FT,N,Nfc} = all([a.corner == b.corner, a.size == b.size])

Base.isless(a::NCube, b::NCube) = Base.isless([a.corner, a.size], [b.corner, b.size])
Base.isequal(a::NCube, b::NCube) = all([a.corner == b.corner, a.size == b.size])
import Base.==
==(a::NCube, b::NCube) = all([a.corner == b.corner, a.size == b.size])

Base.hash(x::NCube, h::UInt) =
    hash((x.corner, x.size), h)

"""
    struct MDBM_Problem{fcT,N,Nf,Nc,Nfc,t01T,t11T,IT,FT}

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
#mutable
struct MDBM_Problem{fcT,N,Nf,Nc,Nfc,t01T,t11T,IT,FT,aT}
    "(memoized) function and constraint in a Tuple"
    fc::fcT
    "The vector of discretized parameters space"
    axes::Axes{N,aT}
    "the bracketing n-cubes (which contains a solution)"
    ncubes::Vector{NCube{IT,FT,N,Nfc}}
    contour_level_fc::Vector #contour levels
    T01::t01T#
    T11pinv::t11T
end

function MDBM_Problem(fc::fcT, axes, ncubes::Vector{<:NCube}, Nf, Nc, Nfc, IT=Int, FT=Float64; contour_level_fc=[nothing, nothing]) where {fcT<:Function}# where {IT<:Integer} where {FT<:AbstractFloat}
    N = length(axes)
    T01 = T01maker(Val(N))
    T11pinv = T11pinvmaker(Val(N))

    println("3 contour_level_fc ", contour_level_fc)
    MDBM_Problem{fcT,N,Nf,Nc,Nfc,typeof(T01),typeof(T11pinv),IT,FT,typeof((axes...,))}(fc, Axes(axes...),
        sort!([NCube{IT,FT,N,Nfc}(MVector{N,IT}([x...]), MVector{N,IT}(ones(IT, length(x))),
            PositionTree(zeros(FT, length(x))), true, MMatrix{N,Nfc,FT}(undef)) for x in Iterators.product((x -> 1:(length(x.ticks)-1)).(axes)...,)][:]), contour_level_fc, T01, T11pinv)
end

function MDBM_Problem(f::Function, axes0::AbstractVector{<:Axis}; constraint::Function=(x...,) -> nothing, memoization::Bool=true,    #Nf=length(f(getindex.(axes0,1)...)),
    Nf=f(getindex.(axes0, 1)...) === nothing ? 0 : length(f(getindex.(axes0, 1)...)),
    Nc=constraint(getindex.(axes0, 1)...) === nothing ? 0 : length(constraint(getindex.(axes0, 1)...)), contour_level_fc=nothing)#Float16(1.), nothing
    Nfc = Nf + Nc
    axes = deepcopy.(axes0)
    argtypesofmyfunc = map(x -> typeof(x).parameters[1], axes)#Argument Type
    AT = Tuple{argtypesofmyfunc...}
    type_f = Base.return_types(f, AT)
    if length(type_f) == 0
        error("input of the function is not compatible with the provided axes")
    else
        RTf = type_f[1]#Return Type of f
    end

    type_con = Base.return_types(constraint, AT)
    if length(type_con) == 0
        error("input of the constraint function is not compatible with the provided axes")
    else
        RTc = type_con[1]#Return Type of the constraint function
    end
    println("Checking memoization ...")
    
    if typeof(f) <: MemF
        println("The function is already memoized: direct useage")
        fun = f
        memoization=true
    else
        println("Creating function with memoization = ", memoization)
        if memoization
            fun = MemF(f, constraint, SortedCache{AT,Tuple{RTf,RTc}}(), RTf, RTc, AT)#Array{MDBMcontainer{RTf,RTc,AT}}(undef, 0))
        else
            #   fun = (x::AT) -> (f(x...)::RTf, constraint(x...)::RTc)::Tuple{RTf,RTc}
            fun = (x) -> (f(x...), constraint(x...))
        end
    end
    Ndim = length(axes)
    if contour_level_fc == nothing
        contour_level_fc = [makezero(RTf), makezero(RTc)]
    end
    MDBM_Problem(fun, axes, Vector{NCube{Int64,Float64,Ndim,Nfc}}(undef, 0), Nf, Nc, Nfc, contour_level_fc=contour_level_fc)
end

function MDBM_Problem(f::Function, a::AbstractVector{<:AbstractVector}; constraint::Function=(x...,) -> nothing, memoization::Bool=true,    #Nf=length(f(getindex.(axes0,1)...)),
    Nf=f(getindex.(a, 1)...) === nothing ? 0 : length(f(getindex.(a, 1)...)),
    Nc=constraint(getindex.(a, 1)...) === nothing ? 0 : length(constraint(getindex.(a, 1)...)), contour_level_fc=nothing)
    axes = [Axis(ax) for ax in a]
    MDBM_Problem(f, axes, constraint=constraint, memoization=memoization, Nf=Nf, Nc=Nc, contour_level_fc=contour_level_fc)#,Vector{NCube{Int64,Float64,Val(Ndim)}}(undef, 0))
end

#notes this keeps the previouse function evaluations
"""
    recreate(mdbm::MDBM_Problem; contour_level_fc=[nothing, nothing], axes_new=mdbm.axes)

Create a new `MDBM_Problem` sharing the memoized function evaluation cache of the input `mdbm`.
This is useful for changing contour levels or axes without re-evaluating the function at existing points.

# Arguments
- `mdbm`: The existing MDBM problem.
- `contour_level_fc`: (Optional) New contour levels.
- `axes_new`: (Optional) New axes definition. Must match the dimension of the original problem.
"""
function recreate(mdbm::MDBM_Problem{fcT,N,Nf,Nc,Nfc,t01T,t11T,IT,FT,aT}; contour_level_fc=[nothing, nothing],axes_new=mdbm.axes) where {fcT,N,Nf,Nc,Nfc,t01T,t11T,IT,FT,aT}
   if N !== length(axes_new) 
    error("The number of axes must be equal to the dimension N of the MDBM problem.\n  N = $(N), length(axes) = $(length(axes_new))")
   end
   axes=Axes(axes_new...)
   mdbm.fc.memoryacc[]=0 #reset the memory acc
    MDBM_Problem{fcT,N,Nf,Nc,Nfc,t01T,t11T,IT,FT,aT}(mdbm.fc, Axes(axes...),
        sort!([NCube{IT,FT,N,Nfc}(MVector{N,IT}([x...]), MVector{N,IT}(ones(IT, length(x))),
            PositionTree(zeros(FT, length(x))), true, MMatrix{N,Nfc,FT}(undef)) for x in Iterators.product((x -> 1:(length(x.ticks)-1)).(axes)...,)][:]), contour_level_fc, mdbm.T01, mdbm.T11pinv)
end

function Base.show(io::IO, mdbm::MDBM_Problem{fcT,N,Nf,Nc,Nfc,t01T,t11T,IT,FT,aT}) where {fcT,N,Nf,Nc,Nfc,t01T,t11T,IT,FT,aT}
    println(io, "Multi-Dimensional Bisection Method Problem")
    println(io, "  parameter dimension:   ", N)#typeof(mdbm).parameters[2])#N
    #println(io, "  co-dimension:          ", Nf)#typeof(mdbm).parameters[3])#Nf
    #println(io, "  number of constraints: ", Nc)#typeof(mdbm).parameters[4])#Nc
    println(io, "  number of co-dimension + constraints:          ", Nfc)#typeof(mdbm).parameters[3])#Nfc
    println(io, "Axes:")
    for k in 1:length(mdbm.axes)
        println(io, "  axis #", k, ": elements: ", length(mdbm.axes[k]), "; elementtype: ", typeof(mdbm.axes[k][1]), "; first: ", mdbm.axes[k][1], "; last: ", mdbm.axes[k][end])
    end

    println(io, "Contour levels: $(mdbm.contour_level_fc)")
    println(io, "number of bracketing n-cubes: ", length(mdbm.ncubes))
    println()

    if (typeof(mdbm.fc) <: MemF)
        println(io, "number of function evaluation: ", length(mdbm.fc.fvalarg))
        println(io, "number of memoized function call: ", mdbm.fc.memoryacc[])
        ratio = length(mdbm.fc.fvalarg) / prod(length.(mdbm.axes))
        if ratio > 0.0
            println(io, "Ratio of function evaluation compared to 'brute-force method': ", ratio)
        end
    else
        println(io, "non-memoized version")
    end
end
