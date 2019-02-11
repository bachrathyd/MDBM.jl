# MDBM__types.jl

using StaticArrays #TODO:(sok helyen lehetne ez, főképp a memoizationban!!!)


struct MDBMcontainer{RTf,RTc,AT}
    funval::RTf
    cval::RTc
    callargs::AT
end


# function (fs::Vector{<:Function})(arg...,)
#     [f(arg...,) for f in fs]
# end
#TODO memoization ofr multiple functions Vector{Function}
struct MemF{RTf,RTc,AT} <:Function
    f::Function
    c::Function
    fvalarg::Vector{MDBMcontainer{RTf,RTc,AT}}#funval,callargs
    memoryacc::Vector{Int64} #TODO: törölni, ha már nem kell
    MemF(f::Function,c::Function,cont::Vector{MDBMcontainer{RTf,RTc,AT}}) where {RTf,RTc,AT}=new{RTf,RTc,AT}(f,c,cont,[Int64(0)])
end

#------------SH version - hidden 'Any'-------------
(memfun::MemF{RTf,RTc,AT})(::Type{RTf},::Type{RTc},args...,) where {RTf,RTc,AT} =( memfun.f(args...,)::RTf, memfun.c(args...,)::RTc)

function (memfun::MemF{RTf,RTc,AT})(args...,) where {RTf,RTc,AT}
    #println("As a normal function call")
    # println(args)
    location=searchsortedfirst(memfun.fvalarg,args,lt=(x,y)->isless(x.callargs,y));
    # println(location)
    # println(length(memfun.fvalarg))
    if length(memfun.fvalarg)<location
        # println("ujraszamolas a végére")
        x=memfun(RTf,RTc,args...,);
        push!(memfun.fvalarg,MDBMcontainer{RTf,RTc,AT}(x...,args))
        return x
    elseif  memfun.fvalarg[location].callargs!=args
        # println("ujraszamolas közé illeszt")
        # println(memfun.fvalarg[location].callargs)
        # println(args)
        x=memfun(RTf,RTc,args...,);
        insert!(memfun.fvalarg, location, MDBMcontainer{RTf,RTc,AT}(x...,args))
        return x
    else
        # println("mar megvolt")
        memfun.memoryacc[1]+=1;
        return (memfun.fvalarg[location].funval,memfun.fvalarg[location].cval);
    end
end
function (memfun::MemF{RTf,RTc,AT})(args::Tuple) where {RTf,RTc,AT}
    #println("Function call with a Tuple---------")
    memfun(args...,)
end



# TODO: Ezt, hogy lehetne meghívni?
function (memfun::MemF{RTf,RTc,AT})(argsVect::Vector{Tuple}) where {RTf,RTc,AT}
    println("MEMO eval for many Tuples")
    [memfun(args) for args in argsVect]
    # c3=sort!(ccc)
    # unique(ccc)
    # cccunique=Set(ccc)
end

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
function fncreator(axes)
    fn="function (ax::$(typeof(axes)))(ind...)\n("
    for i in eachindex(axes)
        fn*="ax[$i][ind[$i]]"
        if i < length(axes)
            fn*=", "
        end
    end
    return fn*")\nend"
end

function createAxesGetindexFunction(axes)
     eval(Meta.parse(fncreator(axes)))
end

function fncreatorTuple(axes)
    fn="function (ax::$(typeof(axes)))(ind)\n("
    for i in eachindex(axes)
        fn*="ax[$i][ind[$i]]"
        if i < length(axes)
            fn*=", "
        end
    end
    return fn*")\nend"
end
function createAxesGetindexFunctionTuple(axes)
     eval(Meta.parse(fncreatorTuple(axes)))
end






# #
# function (axes::Vector{Axis})(is...)
#     getindex.(axes,is)
#     # [mdbm.axes[i][iv] for (i,iv) in enumerate(is)]
# end

#
# function (ax::Vector{Axis})(Coords::Vector{Integer})
#     try
#         # println(i)
#         [ax(i).ticks[Coords[i]] for i in 1:length(ax)]
#         # [ax(i).ticks[Coords[i,ki]] for i in 1:length(ax), ki in 1:size(Coords,2)]
#
#
# # for i in 1:length(ax)
# #     for ki in 1:size(Coords,2)
# #         ax(i).ticks[Coords[i,ki]]
# #     end
# # end
#         # map((x,y)->x.ticks[y], axes,i)
#     catch
#         println("Baj van!!")
#         # println("the number of indexes ($(length(Coords))) in not compatible with the length of axes ($(length(ax)))")
#         # println(Coords)
#         # ((x)->println(length(x.ticks))).(ax)
#     end
# end

function axdoubling!(ax::Axis)
        sort!(append!(ax.ticks, ax.ticks[1:end - 1] + diff(ax.ticks) / 2))
end

# function (axes::Vector{Axis})(i...,)
#     try
#         # println(i)
#         map((x,y)->x.ticks[y], axes,i)
#     catch
#         println("FUCK!")
#         println("the number of indexes ($(length(i))) in not compatible with the length of axes ($(length(axes)))")
#         println(i)
#         ((x)->println(length(x.ticks))).(axes)
#     end
# end

struct NCube{IT<:Integer,FT<:AbstractFloat,ValNdim}
    corner::MVector{ValNdim,IT} #"bottom-left" #Integer index of the axis
    size::MVector{ValNdim,IT}#Integer index of the axis
    posinterp::MVector{ValNdim,FT}#relative coordinate within the cube "(-1:1)" range
    bracketingncube::Bool
    # gradient ::MVector{MVector{T}}
    # curvnorm::Vector{T}
end

function corner(nc::NCube{IT,FT,Ndim},T01)::Vector{MVector} where IT where FT where Ndim
    [nc.corner .+ nc.size .* T for T in T01]
end


function corner(ncubes::Vector{NCube},T01)::Vector{Vector{MVector}}# where IT where FT where Ndim
    [corner(nc,T01) for nc in ncubes]
    # [nc.corner .+ nc.size .* T for nc in ncubes for T in T01]
end


Base.isless(a::NCube{IT,FT,N},b::NCube{IT,FT,N}) where IT where FT where N = Base.isless([a.corner,a.size],[b.corner,b.size])
# Base.isequal(a::NCube{IT,FT,N},b::NCube{IT,FT,N}) where IT where FT where N = all([a.corner==b.corner,a.size==b.size])
import Base.==
==(a::NCube{IT,FT,N},b::NCube{IT,FT,N}) where IT where FT where N = all([a.corner==b.corner,a.size==b.size])


@generated twopow(::Val{n}) where n = 2^n
function T01maker(valk::Val{kdim}) where {kdim}
    # T01=[isodd(x÷(2^y)) for y in 0:(kdim-1), x in 0:(2^kdim-1)]
    # SMatrix{kdim,twopow(valk)}(T01)
        SVector{twopow(valk)}([
        SVector{kdim}([isodd(x÷(2^y)) for y in 0:(kdim-1)])
         for  x in 0:(2^kdim-1)])
end
# //TODO: SH, minek a felső?!?!?!
function T01maker3(valk::Val{kdim}) where {kdim}
    # T01=[isodd(x÷(2^y)) for y in 0:(kdim-1), x in 0:(2^kdim-1)]
    # SMatrix{kdim,twopow(valk)}(T01)
        SVector{kdim^2}([
        SVector{kdim}([isodd(x÷(2^y)) for y in 0:(kdim-1)])
         for  x in 0:(2^kdim-1)])
end

function T11pinvmaker(valk::Val{kdim}) where {kdim}
    T11pinv=([isodd(x÷(2^y)) for y in 0:kdim , x in 0:(2^kdim-1)]*2.0 .-1.0)/(2^kdim)
    SMatrix{kdim+1,twopow(valk)}(T11pinv)
end


struct MDBM_Problem{N}
    fc::Function
    axes::NTuple{N,Axis}
    ncubes::Vector{NCube}
    T01::AbstractVector{<:AbstractVector}#SArray#
    T11pinv::SMatrix

    function MDBM_Problem(fc::Function,axes,ncubes::Vector{NCube{IT,FT,NdimCube}}) where IT<:Integer where FT<:AbstractFloat where NdimCube
        # function MDBM_Problem(fc,axes,ncubes)
        Ndim=length(axes)
        if NdimCube!=Ndim
            error("axis size is not compatible to the nCube size") # internal check, it should never happend TODO: remove this line
        end
        # T01=reshape([isodd(x÷(2^y)) for x in 0:(2^Ndim-1) for y in 0:(Ndim-1)],Ndim,(2^Ndim))
        # T11=reshape([isodd(x÷(2^y)) for x in 0:(2^Ndim-1) for y in 0:Ndim],Ndim+1,(2^Ndim))*2.0 .-1.0
        # T11pinv=T11/(2^Ndim)
        T01=T01maker(Val(Ndim))
        T11pinv=T11pinvmaker(Val(Ndim))
        #new(f,c,axes,[NCube(IT.([x...]),ones(IT,length(x)),zeros(FT,length(axes))) for x in Iterators.product((x->1:(length(x.ticks)-1)).(axes)...,)][:])
        createAxesGetindexFunction((axes...,))
        createAxesGetindexFunctionTuple((axes...,))
        new{Ndim}(fc,(axes...,),
        [NCube{IT,FT,Ndim}(SVector{Ndim,IT}([x...]),SVector{Ndim,IT}(ones(IT,length(x))),SVector{Ndim,FT}(zeros(IT,length(x))),true) for x in Iterators.product((x->1:(length(x.ticks)-1)).(axes)...,)][:]
        ,T01,T11pinv)
    end
end


# function MDBM_Problem{IT,FT}(f::Function, axes::Vector{Axis};constraint::Function=(x...,)->Float16(1.))
#TODO: ha nincs megadva constraint akkor arra ne csináljon memoization (mert biztosan lassabb lesz!)
function MDBM_Problem(f::Function, axes::Vector{<:Axis};constraint::Function=(x...,)->true, memoization::Bool=true)#Float16(1.), nothing
    argtypesofmyfunc=map(x->typeof(x).parameters[1], axes);#Argument Type
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
    fun=MemF(f,constraint,Array{MDBMcontainer{RTf,RTc,AT}}(undef, 0));
    #cons=MemF(constraint,Array{MDBMcontainer{RTc,AT}}(undef, 0));
else
    fun=(x)->(f(x...),constraint(x...));
    #cons=constraint;
end
    #TODO: valami ilyesmi a vektorizációhoz
    Ndim=length(axes)
    MDBM_Problem(fun,axes,Vector{NCube{Int64,Float64,Ndim}}(undef, 0))
end

# function MDBM_Problem(f::Function, a::Vector{<:AbstractVector};constraint::Function=(x...,)->true, memoization::Bool=true)
#     axes=[Axis(ax) for ax in a]
#     argtypesofmyfunc=map(x->typeof(x).parameters[1], axes);#Argument Type
#     AT=Tuple{argtypesofmyfunc...};
#     type_f=Base.return_types(f,AT)
#     if length(type_f)==0
#         error("The input of the function is not compatible with the provided axes")
#     else
#         RTf=type_f[1];#Return Type of f
#     end
#
#     type_con=Base.return_types(constraint,AT)
#     if length(type_con)==0
#         error("The input of the constraint function is not compatible with the provided axes")
#     else
#         RTc=type_con[1];#Return Type of the constraint function
#     end
#     if memoization
#         fun=MemF(f,constraint,Array{MDBMcontainer{RTf,RTc,AT}}(undef, 0));
#         #cons=MemF(constraint,Array{MDBMcontainer{RTc,AT}}(undef, 0));
#     else
#         fun=(x)->(f(x...),constraint(x...));
#         #cons=constraint;
#     end
#     Ndim=length(axes)
#     MDBM_Problem{Ndim}(fun,axes)#,Vector{NCube{Int64,Float64,Val(Ndim)}}(undef, 0))
# end


function corner(mdbm::MDBM_Problem{N})::Vector{Vector{SVector}} where N
    # [corner(nc,mdbm.T01) for nc in mdbm.ncubes]
    [corner(nc,mdbm.T01) for nc in mdbm.ncubes]
    # [nc.corner .+ nc.size .* T for nc in mdbm.ncubes for T in mdbm.T01]
    # c3=sort!(ccc)
    # unique(ccc)
    # cccunique=Set(ccc)
end


function issingchange(FunTupleVector::Vector)::Bool
all(
    [
        [
        any((c)->!isless(c[1][fi],zero(c[1][fi])),FunTupleVector)
        for fi in 1:length(FunTupleVector[1][1])
        ]#any positive (or zeros) value (!isless)

        [
        any((c)->!isless(zero(c[1][fi]),c[1][fi]),FunTupleVector)
        for fi in 1:length(FunTupleVector[1][1])
        ]#anany negative (or zero) value

        [
        any((c)->!isless(c[2][fi],zero(c[2][fi])),FunTupleVector)
        for fi in 1:length(FunTupleVector[1][2])
        ]#any positive (or zeros) value (!isless)
    ]
)#check for all condition
end

#TODO: vagy kockánként kellene csinálni?
function _interpolate!(ncubes::Vector{NCube},mdbm::MDBM_Problem,::Type{Val{0}})
    Ndim=length(mdbm.axes)
    isbracketing=map(nc->issingchange(getcornerval(nc,mdbm)),ncubes)
    deleteat!(ncubes,.!isbracketing)#removeing the non-bracketing ncubes
    # filter!(nc->!issingchange(getcornerval(nc,mdbm)),mdbm.ncubes)

    # ((nc)->nc.posinterp[:].=zero(typeof(nc.posinterp).parameters[1])).(mdbm.ncubes) #set the interpolated relative position to zero
    # ((nc)->nc.posinterp.*=0.0).(mdbm.ncubes) #set the interpolated relative position to zero
    map(nc->nc.posinterp[:].=zeros(typeof(nc.posinterp[1]),Ndim), ncubes)
    return nothing
end


function _interpolate!(ncubes::Vector{NCube},mdbm,::Type{Val{1}})
    Ndim=length(mdbm.axes)
    # Ndim=length(mdbm.axes)
    # # TAntansp=A=hcat([[-1,x...] for x in Iterators.product([(-1.0,1.0) for k in 1:k]...)][:]...);
    # # #     TAn2=inv(TAn.'*TAn);
    # # #     TAtrafo=TAn2*TAn.';
    # # TAtrafo = (TAntansp* transpose(TAntansp) ) \ TAntansp
    # TAtrafoSHORT=hcat([[-1,x...] for x in Iterators.product([(-1.0,1.0) for k in 1:Ndim]...)][:]...)./(2^Ndim);
    #
    for nc in ncubes
        FunTupleVector=getcornerval(nc,mdbm)

        if all([
            any((c)->!isless(c[2][fi],zero(c[2][fi])),FunTupleVector)
            for fi in 1:length(FunTupleVector[1][2])
            ])# do wh have to compute at all?!?!?! ()
            #
            TF=typeof(nc).parameters[2]
            As = Vector{TF}(undef,0)
            ns = Vector{SVector{Ndim,TF}}(undef,0)

            #for f---------------------
            for kf=1:length(FunTupleVector[1][1])
                solloc=mdbm.T11pinv*[FunTupleVector[kcubecorner][1][kf] for kcubecorner=1:length(FunTupleVector)]
                push!(As,solloc[end]);#it is not a real distance within the n-cube (it is ~n*A)!!!
                push!(ns,solloc[1:end-1])
            end

            #for c---if needed: that is- close to the boundary--------
            for kf=1:length(FunTupleVector[1][2])
                #TODO: TODO: mivan a többszörös C teljesülése esetén!?!??!
                if any((c)->!isless(zero(c[2][kf]),c[2][kf]),FunTupleVector) && length(As)<Ndim  # use a constraint till it reduces the dimension to zero (point) and no further
                    solloc=mdbm.T11pinv*[FunTupleVector[kcubecorner][2][kf] for kcubecorner=1:length(FunTupleVector)]
                    push!(As,solloc[end]);#it is not a real distance within the n-cube (it is ~n*A)!!!
                    push!(ns,solloc[1:end-1])
                end
            end

            nsMAT=hcat(ns...)
            #nc.posinterp[:] .= nsMAT * ((transpose(nsMAT) * nsMAT) \ As);
            nc.posinterp[:] .= nsMAT * (inv(transpose(nsMAT) * nsMAT) * As);
            #for c---------------------
        else
            nc.posinterp[:] .=1000.0;#put it outside the cube!
        end
    end

    #TODO: what if it falls outside of the n-cube
    #TODO: it should be removed ->what shall I do with the bracketing cubes?
     #filter!((nc)->norm(nc.posinterp,20.0)>2.0 ,mdbm.ncubes) LinearAlgebre is needed
     #filter!((nc)->sum(nc.posinterp.^20.0)<(2.0 ^20.0),mdbm.ncubes)#1e6~=(2.0 ^20.0)
     filter!((nc)->sum(nc.posinterp.^4.0)<(3.0 ^4.0),mdbm.ncubes)#1e6~=(2.0 ^20.0)
     #filter!((nc)->!any(isnan.(nc.posinterp)),mdbm.ncubes)

    return nothing
end

function _interpolate!(ncubes::Vector{NCube{IT,FT,N}} where IT where FT where N,mdbm,::Type{Val{Ninterp}}) where Ninterp
    error("order $(Ninterp) interpolation is not supperted (yet)")
end

function interpolate!(mdbm::MDBM_Problem;interpolationorder::Int=1)
    _interpolate!(mdbm.ncubes,mdbm, Val{interpolationorder})
end

function getcornerval(mdbm::MDBM_Problem)#get it for all
    #getcornerval.(mdbm.ncubes,Ref(mdbm))
    map((mdbm.fc) ∘ (mdbm.axes),Base.Iterators.flatten(corner(mdbm)))
end

function getcornerval(ncubes::NCube{IT,FT,N} where IT where FT where N,mdbm::MDBM_Problem)
    #PostionsVectors=(mdbm.axes.(corner(nc,mdbm.T01)))
    # (mdbm.fc).( (mdbm.axes).(corner(nc,mdbm.T01)) )
    map((mdbm.fc) ∘ (mdbm.axes),corner(ncubes,mdbm.T01))
end

function doubling!(mdbm::MDBM_Problem,directions::Vector{T}) where T<:Integer
    axdoubling!.(mdbm.axes)
    for nc in mdbm.ncubes
        for dir in directions
            nc.corner[dir]=(nc.corner[dir] -1) *2 +1
            nc.size[dir]*=2
        end
    end
end

function refine!(mdbm::MDBM_Problem{N}; directions::Vector{T}=collect(Int64,1:N)) where N where T<:Integer
    doubling!(mdbm,directions)
    refinencubes!(mdbm.ncubes,directions)
    return nothing
end

function refinencubes!(ncubes::Vector{NCube}, directions::Vector{T}) where T<:Integer #{IT,FT,N} where IT where FT where N
    #TODO:     #if nc.condition>1.0
    # A nem felbontottaknál pedig duplázni kell a felbontást
    for dir in directions
        for nc in ncubes
            nc.size[dir]/=2
        end
        NumofNCubes=length(ncubes)
        append!(ncubes,deepcopy(ncubes))#TODO: ez így kell csinálni?
        for nc in ncubes[1:NumofNCubes]
            nc.corner[dir]=nc.corner[dir]+nc.size[dir]
        end
    end
    sort!(ncubes)#,by=nc->[nc.corner,nc.size])
    return nothing
end


function getinterpolatedpoint(mdbm::MDBM_Problem)
[
    [
        (typeof(mdbm.axes[i].ticks).parameters[1])(
        (mdbm.axes[i].ticks[nc.corner[i]]*(1.0-(nc.posinterp[i]+1.0)/2.0)+
        mdbm.axes[i].ticks[nc.corner[i]+1]*((nc.posinterp[i]+1.0)/2.0))
        )
    for nc in mdbm.ncubes]#::Vector{typeof(mdbm.axes[i].ticks).parameters[1]}
for i in 1:length(mdbm.axes)]
# [
#     [
#         (mdbm.axes[i].ticks[nc.corner[i]]*(typeof(mdbm.axes[i].ticks[1]).(1.0-(nc.posinterp[i]+1.0)/2.0))+
#         mdbm.axes[i].ticks[nc.corner[i]+1]*(typeof(mdbm.axes[i].ticks[1]).((nc.posinterp[i]+1.0)/2.0)))
#     for nc in mdbm.ncubes]::Vector{typeof(mdbm.axes[i].ticks[1])}
# for i in eachindex(mdbm.axes)]
end


function checkneighbour!(mdbm::MDBM_Problem;interpolationorder::Int=1)#only for unite size cubes
allneighind=0:2^length(mdbm.axes)#cornering neightbour also
allneighind=2 .^(0:(length(mdbm.axes)-1))+1#neighbour only in the side
T101=[-mymdbm.T01[allneighind]...,mymdbm.T01[allneighind]...]

compa=(a,b)->all([a.corner==b.corner,a.size==b.size])

ncnew=deepcopy(mdbm.ncubes)
ncubestocheck=deepcopy(mdbm.ncubes)

ismoreneighbour=true
#mTODO: mivan azzal a kockával, aki szomszédként leellenőriztünk, de nem tartalmazott,... majd a követés során újra mellé kerülünk -> így az kétszer lesz ellőnőrizve (ez a Matlabban is rosssz)
while ismoreneighbour
    NumofNCubes=length(ncnew)
    for T in T101

        append!(ncubestocheck,deepcopy(ncnew))#TODO: ez így kell csinálni?
        for nci in 1:NumofNCubes
            ncubestocheck[nci].corner[:]=ncubestocheck[nci].corner+1*T.*ncubestocheck[nci].size
        end
    end
    sort!(ncubestocheck)
    unique(a->[a.corner,a.size], ncubestocheck)
    findin(ncubestocheck,mdbm.ncubes)
end




    sort!(ncubes)#,by=nc->[nc.corner,nc.size])
[
    [
        (typeof(mdbm.axes[i].ticks).parameters[1])(
        (mdbm.axes[i].ticks[nc.corner[i]]*(1.0-(nc.posinterp[i]+1.0)/2.0)+
        mdbm.axes[i].ticks[nc.corner[i]+1]*((nc.posinterp[i]+1.0)/2.0))
        )
    for nc in mdbm.ncubes]#::Vector{typeof(mdbm.axes[i].ticks).parameters[1]}
for i in 1:length(mdbm.axes)]
# [
#     [
#         (mdbm.axes[i].ticks[nc.corner[i]]*(typeof(mdbm.axes[i].ticks[1]).(1.0-(nc.posinterp[i]+1.0)/2.0))+
#         mdbm.axes[i].ticks[nc.corner[i]+1]*(typeof(mdbm.axes[i].ticks[1]).((nc.posinterp[i]+1.0)/2.0)))
#     for nc in mdbm.ncubes]::Vector{typeof(mdbm.axes[i].ticks[1])}
# for i in eachindex(mdbm.axes)]
end


#------------------------



function indexin_sorted(a::Array{T,1}, b::Array{T,1})::Array{Int64,1} where T
        # b must contain all the lements of a
        # a and b must be sorted
    if isempty(a)
        return Array{Int64}(undef,0)
    elseif length(b) == 1 ##much faster!
        return Int64.(a .== b)#*Int64(1)
    else
        out = Array{Int64}(undef, size(a));#zeros(T, size(a))
        leng::Int64 = length(b);

        q::Int64 = 1;
        q1::Int64 = 1;
        q2::Int64 = length(b);
        for k = 1:length(a)

            if (q2 - q1) == 1
                q = (b[q1] == a[k]) ? q1 : q2
            else
                while (b[q] != a[k]) & (q2 > q1 + 1)
                    if b[q] > a[k]
                        q2 = q
                        q =  (q + q1+1) ÷ 2
                    else
                        q1 = q
                        q = (q + q2) ÷ 2
                    end
            #print([q1;q;q2])
            #print([b[q]])
            #println([a[k]])
                end
            end
            out[k] = b[q] == a[k] ? q : 0
            q = max(out[k], q1)
            q1 = q
            q2 = length(b)
          #print("-----")
          #print(out[1:k])
          #println("----")
            if q1 > length(b)#all the element is larger than the last one
                break
            end
        end
        return out
    end
end


# indexin_sorted([-5,0,1,1.1,2,3,4,9.0,11,15,1151],[1.1])
