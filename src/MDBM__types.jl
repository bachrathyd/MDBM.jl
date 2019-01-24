# MDBM__types.jl

using StaticArrays #TODO:(sok helyen lehetne ez, főképp a memoizationban!!!)


struct MDBMcontainer{RT,AT}
    funval::RT
    callargs::AT
end



# function (fs::Vector{<:Function})(arg...,)
#     [f(arg...,) for f in fs]
# end
#TODO memoization ofr multiple functions Vector{Function}
struct MemF{RT,AT} <:Function
    f::Function
    fvalarg::Vector{MDBMcontainer{RT,AT}}#funval,callargs
    memoryacc::Vector{Int64} #TODO: törölni, ha már nem kell
    MemF(f::Function,a::Vector{MDBMcontainer{RT,AT}}) where {RT,AT}=new{RT,AT}(f,a,[Int64(0)])
end


# function MemF{RT,AT}(rT::RT,aT::AT)
#     MemF{RT,AT}(rT,aT,Int64(1))
# end

#------------SH version - hidden 'Any'-------------
(memfun::MemF{RT,AT})(::Type{RT},args...,) where {RT,AT} = memfun.f(args...,)::RT

function (memfun::MemF{RT,AT})(args...,) where {RT,AT}
    location=searchsortedfirst(memfun.fvalarg,args,lt=(x,y)->isless(x.callargs,y));
    if length(memfun.fvalarg)<location
        #println("ujraszamolas a végére")
        x=memfun(RT,args...,);
        push!(memfun.fvalarg,MDBMcontainer{RT,AT}(x,args))
        return x
    elseif  memfun.fvalarg[location].callargs!=args
        #println("ujraszamolas közé illeszt")
        x=memfun(RT,args...,);
        insert!(memfun.fvalarg, location, MDBMcontainer{RT,AT}(x,args))
        return x
    else
        #println("mar megvolt")
        memfun.memoryacc[1]+=1;
        return memfun.fvalarg[location].funval;
    end
end
#-------------------------

# function (memfun::MemF{RT,AT})(Vector{T}) where {T,RT,AT}
#     locations=indexin_sorted(memfun.fvalarg,args,lt=(x,y)->isless(x.callargs,y));
#     for x in locations if x!=0
#         x=memfun(RT,args...,)::RT;
#     end
# end
#
# function indexin_sorted(a::Array{T,1}, b::Array{T,1})::Array{Int64,1} where T
#         # b must contain all the lements of a
#         # a and b must be sorted
#     if isempty(a)
#         return Array{Int64}(undef,0)
#     elseif length(b) == 1 ##much faster!
#         return Int64.(a .== b)#*Int64(1)
#     else
#         out = Array{Int64}(undef, size(a));#zeros(T, size(a))
#         leng::Int64 = length(b);
#
#         q::Int64 = 1;
#         q1::Int64 = 1;
#         q2::Int64 = length(b);
#         for k = 1:length(a)
#
#             if (q2 - q1) == 1
#                 q = (b[q1] == a[k]) ? q1 : q2
#             else
#                 while (b[q] != a[k]) & (q2 > q1 + 1)
#                     if b[q] > a[k]
#                         q2 = q
#                         q =  (q + q1+1) ÷ 2
#                     else
#                         q1 = q
#                         q = (q + q2) ÷ 2
#                     end
#             #print([q1;q;q2])
#             #print([b[q]])
#             #println([a[k]])
#                 end
#             end
#             out[k] = b[q] == a[k] ? q : 0
#             q = max(out[k], q1)
#             q1 = q
#             q2 = length(b)
#           #print("-----")
#           #print(out[1:k])
#           #println("----")
#             if q1 > length(b)#all the element is larger than the last one
#                 break
#             end
#         end
#         return out
#     end
# end
#
#
# # indexin_sorted([-5,0,1,1.1,2,3,4,9.0,11,15,1151],[1.1])




struct Axis{T}
    ticks::Vector{T}
    name
end

Base.getindex(ax::Axis, ind) = ax.ticks[ind]
Base.setindex!(ax::Axis, X, inds...,) = setindex!(ax.ticks, X, inds...,)
Base.size(ax::Axis) = size(ax.ticks)


function Axis(a::Vector{T}) where T
    Axis(a, :unknown)
end

function axdoubling!(ax::Axis)
        sort!(append!(ax.ticks, ax.ticks[1:end - 1] + diff(ax.ticks) / 2))
end


#TODO SH: jó ez így ??
function Axis(a)
    Axis([a...], :unknown)
end
function Axis(a,b)
    Axis([a...], b)
end

# function Axis(a::AbstractVector, name::Symbol=:unknown;dtype::Type=Float64)
#     Axis(convert(Vector{dtype}, a), name)
# end

function Axis(a::Axis)
    a
end




struct NCube{IT<:Integer,FT<:AbstractFloat}
    corner::Vector{IT} #"bottom-left" #Integer index of the axis
    size::Vector{IT}#Integer index of the axis
    posinterp::Vector{FT}#relative coordinate within the cube "(-1:1)" range
    #gradient ::Vector{T}
    # curvnorm::Vector{T}
end

# TOOD: egyszerűsíteni a létrehozást (bár lehet, hogy nem kell, csak 1-2 helyen fog szerepelni (max))
# function NCube{IT<:Integer,FT<:AbstractFloat}(x::Vector{IT}) where IT where FT
#     NCube(IT.(x),ones(IT,length(x)),Vector{FT}(undef,length(x)))
# end




# struct MDBM_Problem
#     f::MemF
#     c::MemF
#     axes::Vector{Axis}
#     ncubes::Vector{NCube}
#     function MDBM_Problem(f,c,axes,ncubes::Vector{NCube{IT,FT}}) where IT where FT
#         #new(f,c,axes,[NCube(IT.([x...]),ones(IT,length(x)),zeros(FT,length(axes))) for x in Iterators.product((x->1:(length(x.ticks)-1)).(mdbmaxes)...,)][:])
#         new(f,c,axes,[NCube(IT.([x...]),ones(IT,length(x)),Vector{FT}(undef,length(axes))) for x in Iterators.product((x->1:(length(x.ticks)-1)).(mdbmaxes)...,)][:])
#     end
# end

struct MDBM_Problem
    f::Function
    c::Function
    axes::Vector{Axis}
    ncubes::Vector{NCube}
    function MDBM_Problem(f,c,axes,ncubes::Vector{NCube{IT,FT}}) where IT where FT
        #new(f,c,axes,[NCube(IT.([x...]),ones(IT,length(x)),zeros(FT,length(axes))) for x in Iterators.product((x->1:(length(x.ticks)-1)).(mdbmaxes)...,)][:])
        new(f,c,axes,[NCube(IT.([x...]),ones(IT,length(x)),Vector{FT}(undef,length(axes))) for x in Iterators.product((x->1:(length(x.ticks)-1)).(mdbmaxes)...,)][:])
    end
end
function (axes::Vector{Axis})(i::Vector)
    map((x,y)->x.ticks[y], mdbm.axes,i)
end

function (axes::Vector{Axis})(i...,)
    try
        map((x,y)->x.ticks[y], mdbm.axes,i)
    catch
        println("the number of indexes ($(length(i))) in not compatible with the length of axes ($(length(axes)))")
        println(i)
        ((x)->println(length(x.ticks))).(mdbm.axes)
    end
end


# function MDBM_Problem{IT,FT}(f::Function, mdbmaxes::Vector{Axis};constraint::Function=(x...,)->Float16(1.))
#TODO: ha nincs megadva constraint akkor arra ne csináljon memoization (mert biztosan lassabb lesz!)
function MDBM_Problem(f::Function, mdbmaxes::Vector{<:Axis};constraint::Function=(x...,)->true, memoization::Bool=true)#Float16(1.)
    argtypesofmyfunc=map(x->typeof(x).parameters[1], mdbmaxes);#Argument Type
    AT=Tuple{argtypesofmyfunc...};
    type_f=Base.return_types(f,AT)
    if length(type_f)==0
        error("input of the function is not compatible with the provided axes")
    else
        RTf=type_f[1];#Return Type of f
    end

    type_con=Base.return_types(constraint,AT)
    if length(type_con)==0
        error("input of the constraint function is not compatible with the provided axes")
    else
        RTc=type_con[1];#Return Type of the constraint function
    end

if memoization
    fun=MemF(f,Array{MDBMcontainer{RTf,AT}}(undef, 0));
    cons=MemF(constraint,Array{MDBMcontainer{RTc,AT}}(undef, 0));
else
    fun=f;
    cons=constraint;
end
    #TODO: valami ilyesmi a vektorizációhoz
    #fun=MemF((x)->[f(x...,),constraint(x...,)],Array{MDBMcontainer{RTc,AT}}(undef, 0));


    MDBM_Problem(fun,cons,mdbmaxes,Vector{NCube{Int64,Float64}}(undef, 0))
end

# {AbstractArray{T,1} where T}
# function MDBM_Problem(f::Function, a::AbstractVector;constraint::Function=(x...,)->Float16(1.))
function MDBM_Problem(f::Function, a::Vector{<:AbstractVector};constraint::Function=(x...,)->true, memoization::Bool=true)
    mdbmaxes=[Axis(ax) for ax in a]
    argtypesofmyfunc=map(x->typeof(x).parameters[1], mdbmaxes);#Argument Type
    AT=Tuple{argtypesofmyfunc...};
    type_f=Base.return_types(f,AT)
    if length(type_f)==0Ű
        error("The input of the function is not compatible with the provided axes")
    else
        RTf=type_f[1];#Return Type of f
    end

    type_con=Base.return_types(constraint,AT)
    if length(type_con)==0
        error("The input of the constraint function is not compatible with the provided axes")
    else
        RTc=type_con[1];#Return Type of the constraint function
    end
    if memoization
        fun=MemF(f,Array{MDBMcontainer{RTf,AT}}(undef, 0));
        cons=MemF(constraint,Array{MDBMcontainer{RTc,AT}}(undef, 0));
    else
        fun=f;
        cons=constraint;
    end
    MDBM_Problem(fun,cons,mdbmaxes,Vector{NCube{Int64,Float64}}(undef, 0))
end




#TODO: vagy kockánként kellene csinálni?
function _interpolate!(mdbm,::Type{Val{0}})
    print("0 cuccc happening")
    issignchangeforeveryfunction=(ncube)->
    all(

        [
            [
            any((c)->isless(c[fi],zero(c[fi])),ncube[1])
            for fi in 1:length(ncube[1][1])
            ]#any positive value

            [
            any((c)->isless(zero(c[fi]),c[fi]),ncube[1])
            for fi in 1:length(ncube[1][1])
            ]#anany positive valueny positive value

            [
            any((c)->isless(zero(c[fi]),c[fi]),ncube[2])
            for fi in 1:length(ncube[2][1])
            ]#if any positive valueny positive value
        ][:]
    )#chack for all
    isbracketing=[issignchangeforeveryfunction(c) for c in getcornerval(mdbm)]
    deleteat!(mdbm.ncubes,.!isbracketing)
    ((nc)->nc.posinterp[:].=zero(typeof(nc.posinterp).parameters[1])).(mdbm.ncubes)
    return nothing
end


function _interpolate!(mdbm,::Type{Val{1}})
    print("1 cuccc happening")
    Ndim=length(mdbm.axes)
    # TAntansp=A=hcat([[-1,x...] for x in Iterators.product([(-1.0,1.0) for k in 1:k]...)][:]...);
    # #     TAn2=inv(TAn.'*TAn);
    # #     TAtrafo=TAn2*TAn.';
    # TAtrafo = (TAntansp* transpose(TAntansp) ) \ TAntansp
    TAtrafoSHORT=A=hcat([[-1,x...] for x in Iterators.product([(-1.0,1.0) for k in 1:Ndim]...)][:]...)./(2^Ndim);
    # norm(TAtrafo-TAtrafoSHORT)



    for nc in mdbm.ncubes
        fcvals=getcornerval(nc,mdbm)

        if all([
            any((c)->isless(zero(c[fi]),c[fi]),fcvals[2])
            for fi in 1:length(fcvals[2][1])
            ])# do wh have to compute at all?!?!?! ()

            TF=typeof(nc).parameters[2]
            As = Vector{TF}(undef,0)
            ns = Vector{SVector{length(mdbm.axes),TF}}(undef,0)

            #for f---------------------
            for kf=1:length(fcvals[1][1])
                solloc=TAtrafoSHORT*[fcvals[1][kcube][kf] for kcube=1:length(fcvals[1])]
                push!(As,solloc[1]);#it is not a real distance within the n-cube (it is ~n*A)!!!
                push!(ns,solloc[2:end]);
            end

            #for c---if needed: that is- close to the boundary--------
            for kf=1:length(fcvals[2][1])
                if any((c)->isless(c[kf],zero(c[kf])),fcvals[2])
                    solloc =TAtrafoSHORT*[fcvals[2][kcube][kf] for kcube=1:length(fcvals[2])]
                    push!(As,solloc[1]);#it is not a real distance within the n-cube (it is ~n*A)!!!
                    push!(ns,solloc[2:end]);
                end
            end
            nsMAT=hcat(ns...)
            #nc.posinterp[:] .= nsMAT * ((transpose(nsMAT) * nsMAT) \ As);
            nc.posinterp[:] .= nsMAT * (inv(transpose(nsMAT) * nsMAT) * As);
            #for c---------------------
        else
            nc.posinterp[:] .=1000.0;
        end
    end

    #TODO: what if it falls outside of the n-cube
    #TODO: it should be removed ->what shall I do with the bracketing cubes?
     #filter!((nc)->norm(nc.posinterp,20.0)>2.0 ,mdbm.ncubes) LinearAlgebre is needed
     filter!((nc)->sum(nc.posinterp.^20.0)<(2.0 ^20.0),mdbm.ncubes)#1e6~=(2.0 ^20.0)

    return nothing
end

function _interpolate!(mdbm,::Type{Val{N}}) where N
    error("order $(N) interpolation is not supperted (yet)")
end

function interpolate!(mdbm::MDBM_Problem;interpolationorder::Int=1)
    _interpolate!(mdbm, Val{interpolationorder})
end

function getcornerval(mdbm::MDBM_Problem)#get it for all
    [
    getcornerval(nc,mdbm)
            for nc in mdbm.ncubes]
end
# function mdbmcall(mdbm::MDBM_Problem,x)
#     [mdbm.f(x...,),mdbm.c(x...,)]
# end

function getcornerval(nc::NCube{IT,FT} where IT where FT,mdbm::MDBM_Problem)
 [
    [
    mdbm.f.((mdbm.axes).(x...,)...,)
    for x in
        Iterators.product(
        [(nc.corner[n],nc.corner[n]+nc.size[n])
        for n in 1:length(nc.corner)]
            ...,)
            ][:]
,
    [
    mdbm.c.((mdbm.axes).(x...,)...,)
    for x in
        Iterators.product(
        [(nc.corner[n],nc.corner[n]+nc.size[n])
        for n in 1:length(nc.corner)]
            ...,)
            ][:]
]
end


function refine!(mdbm::MDBM_Problem)
    axdoubling!.(mdbm.axes)
    [nc.corner[:].=((nc.corner .-1.0) .*2.0 .+1.0 )   for nc in mdbm.ncubes]
    [nc.posinterp[:].=Vector{typeof(nc.posinterp[1])}(undef,length(mdbm.axes))    for nc in mdbm.ncubes]

    append!(mdbm.ncubes,vcat([

       [

       # push!(mdbm.ncubes,NCube(typeof(nc.corner[1]).([x...]),nc.size,Vector{typeof(nc.posinterp[1])}(undef,length(mdbm.axes))))
        #
       NCube(typeof(nc.corner[1]).([x...]),nc.size,Vector{typeof(nc.posinterp[1])}(undef,length(mdbm.axes)))
       for x in
           Iterators.product(
           [(nc.corner[n],nc.corner[n]+nc.size[n])
           for n in 1:length(nc.corner)]
            ...,)
            if [x...]!=nc.corner
        ]
    for nc in mdbm.ncubes]...))
            #TODO:     #if nc.condition>1.0
                 # A nem felbontottaknál pedig duplázni kell a felbontást
return nothing
end


function getinterpolatedpoint(mdbm::MDBM_Problem)
[
    [
        mdbm.axes[i].ticks[nc.corner[i]]*(1-(nc.posinterp[i]+1)/2)+
        mdbm.axes[i].ticks[nc.corner[i]+1]*((nc.posinterp[i]+1)/2)
    for nc in mdbm.ncubes]::Vector{typeof(mdbm.axes[1].ticks[1])}
for i in 1:length(mdbm.axes)]
end
