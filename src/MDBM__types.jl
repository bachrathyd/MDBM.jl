
using StaticArrays #TODO:(sok helyen lehetne ez, főképp a memoizationban!!!)


struct MDBMcontainer{RTf,RTc,AT}
    funval::RTf
    cval::RTc
    callargs::AT
end

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
function (memfun::MemF{RTf,RTc,AT})(args::Tuple) where {RTf,RTc,AT}
    memfun(args...,)
end


# TODO: Ezt, hogy lehetne meghívni? a csoportos futtatáshoz???
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

function axdoubling!(ax::Axis)
        sort!(append!(ax.ticks, ax.ticks[1:end - 1] + diff(ax.ticks) / 2); alg=QuickSort)
end

struct NCube{IT<:Integer,FT<:AbstractFloat,ValNdim}
    corner::MVector{ValNdim,IT} #"bottom-left" #Integer index of the axis
    size::MVector{ValNdim,IT}#Integer index of the axis
    posinterp::MVector{ValNdim,FT}#relative coordinate within the cube "(-1:1)" range
    bracketingncube::Bool
    # gradient ::MVector{MVector{T}}
    # curvnorm::Vector{T}
end

function corner(nc::NCube{IT,FT,N},T01)::Vector{MVector} where IT where FT where N
    [nc.corner .+ nc.size .* T for T in T01]
end


function corner(ncubes::Vector{NCube{IT,FT,N}},T01)::Vector{Vector{MVector}} where IT where FT where N
    [corner(nc,T01) for nc in ncubes]
    # [nc.corner .+ nc.size .* T for nc in ncubes for T in T01]
end


Base.isless(a::NCube{IT,FT,N},b::NCube{IT,FT,N}) where IT where FT where N = Base.isless([a.corner,a.size],[b.corner,b.size])
Base.isequal(a::NCube{IT,FT,N},b::NCube{IT,FT,N}) where IT where FT where N = all([a.corner==b.corner,a.size==b.size])
import Base.==
==(a::NCube{IT,FT,N},b::NCube{IT,FT,N}) where IT where FT where N = all([a.corner==b.corner,a.size==b.size])
==(a::NCube,b::NCube)= all([a.corner==b.corner,a.size==b.size])
==(a::NCube{IT,FT,N},b::NCube) where IT where FT where N = all([a.corner==b.corner,a.size==b.size])
==(a::NCube,b::NCube{IT,FT,N}) where IT where FT where N = all([a.corner==b.corner,a.size==b.size])


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

struct MDBM_Problem{N,Nf,Nc} #where N where Nf where Nc
    fc::Function
    axes::NTuple{N,Axis}
    ncubes::Vector{NCube{IT,FT,N}} where IT<:Integer where FT<:AbstractFloat
    T01::AbstractVector{<:AbstractVector}#SArray#
    T11pinv::SMatrix
end

function MDBM_Problem(fc::Function,axes,ncubes::Vector{NCube{IT,FT,N}},Nf,Nc) where IT<:Integer where FT<:AbstractFloat where N
    # function MDBM_Problem(fc,axes,ncubes)
    T01=T01maker(Val(N))
    T11pinv=T11pinvmaker(Val(N))
    createAxesGetindexFunction((axes...,))
    createAxesGetindexFunctionTuple((axes...,))
    MDBM_Problem{N,Nf,Nc}(fc,(axes...,),
    [NCube{IT,FT,N}(SVector{N,IT}([x...]),SVector{N,IT}(ones(IT,length(x))),SVector{N,FT}(zeros(IT,length(x))),true) for x in Iterators.product((x->1:(length(x.ticks)-1)).(axes)...,)][:]
    ,T01,T11pinv)
end

# function MDBM_Problem{IT,FT}(f::Function, axes::Vector{Axis};constraint::Function=(x...,)->Float16(1.))
#TODO: ha nincs megadva constraint akkor arra ne csináljon memoization (mert biztosan lassabb lesz!)
function MDBM_Problem(f::Function, axes::Vector{<:Axis};constraint::Function=(x...,)->true, memoization::Bool=true,
    Nf=length(f(getindex.(axes,1)...)),
    Nc=length(constraint(getindex.(axes,1)...)))#Float16(1.), nothing

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
    MDBM_Problem(fun,axes,Vector{NCube{Int64,Float64,Ndim}}(undef, 0),Nf,Nc)
end

function MDBM_Problem(f::Function, a::Vector{<:AbstractVector};constraint::Function=(x...,)->true, memoization::Bool=true,
    Nf=length(f(getindex.(a,1)...)),
    Nc=length(constraint(getindex.(a,1)...)))

    axes=[Axis(ax) for ax in a]
    MDBM_Problem(f,axes,constraint=constraint,memoization=memoization,Nf=Nf,Nc=Nc)#,Vector{NCube{Int64,Float64,Val(Ndim)}}(undef, 0))
end

function corner(mdbm::MDBM_Problem{N,Nf,Nc})::Vector{Vector{SVector}} where N where Nf where Nc
    [corner(nc,mdbm.T01) for nc in mdbm.ncubes]
end

function issingchange(FunTupleVector::Vector,Nf::Integer,Nc::Integer)::Bool
all(
    [
        [
        any((c)->!isless(c[1][fi],zero(c[1][fi])),FunTupleVector)
        for fi in 1:Nf#length(FunTupleVector[1][1])
        ]#any positive (or zeros) value (!isless)

        [
        any((c)->!isless(zero(c[1][fi]),c[1][fi]),FunTupleVector)
        for fi in 1:Nf#length(FunTupleVector[1][1])
        ]#anany negative (or zero) value

        [
        any((c)->!isless(c[2][fi],zero(c[2][fi])),FunTupleVector)
        for fi in 1:Nc#length(FunTupleVector[1][2])
        ]#any positive (or zeros) value (!isless)
    ]
)#check for all condition
end

#TODO: vagy kockánként kellene csinálni?
# function _interpolate!(ncubes::Vector{NCube},mdbm::MDBM_Problem{N},::Type{Val{0}})
#     Ndim=length(mdbm.axes)
#     isbracketing=map(nc->issingchange(getcornerval(nc,mdbm)),ncubes)
#     deleteat!(ncubes,.!isbracketing)#removeing the non-bracketing ncubes
#     # filter!(nc->!issingchange(getcornerval(nc,mdbm)),mdbm.ncubes)
#
#     # ((nc)->nc.posinterp[:].=zero(typeof(nc.posinterp).parameters[1])).(mdbm.ncubes) #set the interpolated relative position to zero
#     # ((nc)->nc.posinterp.*=0.0).(mdbm.ncubes) #set the interpolated relative position to zero
#     map(nc->nc.posinterp[:].=zeros(typeof(nc.posinterp[1]),Ndim), ncubes) #TODO: vagy egy forciklussal?
#     return nothing
# end
function _interpolate!(ncubes::Vector{NCube{IT,FT,N}},mdbm::MDBM_Problem{N,Nf,Nc},::Type{Val{0}}) where IT where FT where N where Nf where Nc
    isbracketing=map(nc->issingchange(getcornerval(nc,mdbm),Nf,Nc),ncubes)
    deleteat!(ncubes,.!isbracketing)#removeing the non-bracketing ncubes
    # filter!(nc->!issingchange(getcornerval(nc,mdbm)),mdbm.ncubes)

    # ((nc)->nc.posinterp[:].=zero(typeof(nc.posinterp).parameters[1])).(mdbm.ncubes) #set the interpolated relative position to zero
    # ((nc)->nc.posinterp.*=0.0).(mdbm.ncubes) #set the interpolated relative position to zero
    map(nc->nc.posinterp[:].=zeros(typeof(nc.posinterp[1]),N), ncubes) #TODO: vagy egy forciklussal?
    return nothing
end

function _interpolate!(ncubes::Vector{NCube{IT,FT,N}},mdbm::MDBM_Problem{N,Nf,Nc},::Type{Val{1}}) where IT where FT where N where Nf where Nc
#function _interpolate!(ncubes::Vector{NCube{IT,FT,N}} where IT where FT where N,mdbm,::Type{Val{1}})
    # N=length(mdbm.axes)
    # N=length(mdbm.axes)
    # # TAntansp=A=hcat([[-1,x...] for x in Iterators.product([(-1.0,1.0) for k in 1:k]...)][:]...);
    # # #     TAn2=inv(TAn.'*TAn);
    # # #     TAtrafo=TAn2*TAn.';
    # # TAtrafo = (TAntansp* transpose(TAntansp) ) \ TAntansp
    # TAtrafoSHORT=hcat([[-1,x...] for x in Iterators.product([(-1.0,1.0) for k in 1:N]...)][:]...)./(2^N);
    #

    for nc in ncubes
        FunTupleVector=getcornerval(nc,mdbm)

        if all([
            any((c)->!isless(c[2][fi],zero(c[2][fi])),FunTupleVector)
            for fi in 1:length(FunTupleVector[1][2])
            ])# do wh have to compute at all?!?!?! ()
            #
            # As = Vector{FT}(undef,0)
            # ns = Vector{SVector{N,FT}}(undef,0)
            As = zero(MVector{Nf+Nc,FT})
            ns = zero(MMatrix{N,Nf+Nc,FT})


            #for f---------------------
            for kf=1:Nf#length(FunTupleVector[1][1])
                solloc=mdbm.T11pinv*[FunTupleVector[kcubecorner][1][kf] for kcubecorner=1:length(FunTupleVector)]
                As[kf]=solloc[end];
                ns[:,kf].=solloc[1:end-1];
                # push!(As,solloc[end]);#it is not a real distance within the n-cube (it is ~n*A)!!!
                # push!(ns,solloc[1:end-1])
            end

            #for c---if needed: that is- close to the boundary--------

            activeCostraint=0;
            for kf=1:Nc#length(FunTupleVector[1][2])
                #TODO: TODO: mivan a többszörös C teljesülése esetén!?!??!
                if any((c)->!isless(zero(c[2][kf]),c[2][kf]),FunTupleVector) && length(As)<N  # use a constraint till it reduces the dimension to zero (point) and no further
                    solloc=mdbm.T11pinv*[FunTupleVector[kcubecorner][2][kf] for kcubecorner=1:length(FunTupleVector)]
                    activeCostraint+=1;
                    As[Nf+activeCostraint]=solloc[end];
                    ns[Nf+activeCostraint,kf].=solloc[1:end-1];
                    # push!(As,solloc[end]);#it is not a real distance within the n-cube (it is ~n*A)!!!
                    # push!(ns,solloc[1:end-1]) nc.posinterp[:] .= ns
                end
            end

            if (Nf+activeCostraint)==0
                nc.posinterp[:] .=0.0
            else
                # nsMAT=hcat(ns...)
                # # nc.posinterp[:] .= nsMAT * ((transpose(nsMAT) * nsMAT) \ As);
                # nc.posinterp[:] .= nsMAT * (inv(transpose(nsMAT) * nsMAT) * As);
                # # nc.posinterp[:] .= transpose(nsMAT) \ As;
                nsMAT=ns[:,1:(Nf+activeCostraint)];
                nc.posinterp[:] .= nsMAT * (inv(transpose(nsMAT) * nsMAT) * As[1:(Nf+activeCostraint)]);
                # nc.posinterp[:] .= ns[:,1:(Nf+activeCostraint)] \As[1:(Nf+activeCostraint)];#TEST
            end

            #for c---------------------
        else
            nc.posinterp[:] .=1000.0;#put it outside the cube!
        end
    end

    #TODO: what if it falls outside of the n-cube
    #TODO: it should be removed ->what shall I do with the bracketing cubes?
     #filter!((nc)->norm(nc.posinterp,20.0)>2.0 ,mdbm.ncubes) LinearAlgebre is needed
     #filter!((nc)->sum(nc.posinterp.^20.0)<(2.0 ^20.0),mdbm.ncubes)#1e6~=(2.0 ^20.0)
     filter!((nc)->sum((abs.(nc.posinterp)).^10.0)<(1.2 ^10.0),mdbm.ncubes)#1e6~=(2.0 ^20.0)
     #filter!((nc)->!any(isnan.(nc.posinterp)),mdbm.ncubes)

    return nothing
end

function _interpolate!(ncubes::Vector{NCube{IT,FT,N}},mdbm::MDBM_Problem{N,Nf,Nc},::Type{Val{Ninterp}}) where IT where FT where N where Ninterp where Nf where Nc
    error("order $(Ninterp) interpolation is not supperted (yet)")
end

function interpolate!(mdbm::MDBM_Problem{N,Nf,Nc};interpolationorder::Int=1) where N where Nf where Nc
    _interpolate!(mdbm.ncubes,mdbm, Val{interpolationorder})
end


function axesextend!(mdbm::MDBM_Problem{N,Nf,Nc},axisnumber::Integer;prepend::AbstractArray=[],append::AbstractArray=[]) where N where Nf where Nc
    preplength=length(prepend);
    if preplength>0
        prepend!(mdbm.axes[axisnumber].ticks,prepend)
        for nc in mdbm.ncubes
            nc.corner[axisnumber]+=preplength;
        end
    end
    if length(append)>0
        append!(mdbm.axes[axisnumber].ticks,append)
    end
end

function getcornerval(mdbm::MDBM_Problem{N,Nf,Nc}) where N where Nf where Nc#get it for all
    #getcornerval.(mdbm.ncubes,Ref(mdbm))
    map((mdbm.fc) ∘ (mdbm.axes),Base.Iterators.flatten(corner(mdbm)))
end

function getcornerval(ncubes::NCube{IT,FT,N},mdbm::MDBM_Problem{N,Nf,Nc}) where IT where FT where N where Nf where Nc
    #PostionsVectors=(mdbm.axes.(corner(nc,mdbm.T01)))
    # (mdbm.fc).( (mdbm.axes).(corner(nc,mdbm.T01)) )
    map((mdbm.fc) ∘ (mdbm.axes),corner(ncubes,mdbm.T01))
end

function doubling!(mdbm::MDBM_Problem{N,Nf,Nc},directions::Vector{T}) where T<:Integer  where N where Nf where Nc
    axdoubling!.(mdbm.axes[directions])
    for nc in mdbm.ncubes
        for dir in directions
            nc.corner[dir]=(nc.corner[dir] -1) *2 +1
            nc.size[dir]*=2
        end
    end
end

function refine!(mdbm::MDBM_Problem{N,Nf,Nc}; directions::Vector{T}=collect(Int64,1:N)) where N where Nf where Nc where T<:Integer
    doubling!(mdbm,directions)
    refinencubes!(mdbm.ncubes,directions)
    return nothing
end

function refinencubes!(ncubes::Vector{NCube{IT,FT,N}}, directions::Vector{T}) where IT where FT where N where T<:Integer #{IT,FT,N} where IT where FT where N
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
    sort!(ncubes; alg=QuickSort)#,by=nc->[nc.corner,nc.size])
    return nothing
end


function getinterpolatedpoint(mdbm::MDBM_Problem{N,Nf,Nc}) where N where Nf where Nc
[
    [
        (typeof(mdbm.axes[i].ticks).parameters[1])(
        (mdbm.axes[i].ticks[nc.corner[i]]*(1.0-(nc.posinterp[i]+1.0)/2.0)+
        mdbm.axes[i].ticks[nc.corner[i]+1]*((nc.posinterp[i]+1.0)/2.0))
        )
    for nc in mdbm.ncubes]#::Vector{typeof(mdbm.axes[i].ticks).parameters[1]}
for i in 1:length(mdbm.axes)]
end
function getinterpolatedpoint(nc::NCube{IT,FT,N},mdbm::MDBM_Problem{N,Nf,Nc}) where IT where FT where N where Nf where Nc
[
        (typeof(mdbm.axes[i].ticks).parameters[1])(
        (mdbm.axes[i].ticks[nc.corner[i]]*(1.0-(nc.posinterp[i]+1.0)/2.0)+
        mdbm.axes[i].ticks[nc.corner[i]+1]*((nc.posinterp[i]+1.0)/2.0))
        )#::Vector{typeof(mdbm.axes[i].ticks).parameters[1]}
for i in 1:length(mdbm.axes)]
end

# function generateneighbours(ncubes::Vector{NCube},mdbm::MDBM_Problem{N,Nf,Nc}) where N  where Nf where Nc
function generateneighbours(ncubes::Vector{NCube{IT,FT,N}},mdbm::MDBM_Problem{N,Nf,Nc}) where IT where FT where N  where Nf where Nc
    #-------all faceing/cornering neightbours----------------
    # #neighbourind=1:2^length(mdbm.axes) #cornering neightbours also - unnecessary
    # neighbourind=2 .^(0:(length(mdbm.axes)-1)) .+ 1 #neighbour only on the side
    # T101=[-mdbm.T01[neighbourind]...,mdbm.T01[neighbourind]...]
    #
    # nc_neighbour = Array{typeof(mdbm.ncubes[1])}(undef,0)
    # NumofNCubes=length(ncubes)
    # for iT in 1:length(T101)
    #     append!(nc_neighbour,deepcopy(ncubes))#TODO: ez így kell csinálni?
    #     for nci in ((1+NumofNCubes*(iT-1)):(NumofNCubes+NumofNCubes*(iT-1)))
    #         nc_neighbour[nci].corner[:]=nc_neighbour[nci].corner+T101[iT].*nc_neighbour[nci].size
    #     end
    # end
    #-------all faceing/cornering neightbours----------------


    #-------direcational neightbours----------------
    nc_neighbour = Array{typeof(mdbm.ncubes[1])}(undef,0)
    Nface=2^N
    indpos=[[mdbm.T01[i][dir] for i in 1:Nface] for dir in 1:N]
    indneg=[[!mdbm.T01[i][dir] for i in 1:Nface] for dir in 1:N]

    for nc in ncubes
        fcvals=getcornerval(nc,mdbm)
        for dir in 1:N
            # indpos=[mdbm.T01[i][dir] for i in 1:Nface]
            # indneg=[!mdbm.T01[i][dir] for i in 1:Nface]
            if issingchange(fcvals[indpos[dir]],Nf,Nc)
                push!(nc_neighbour,deepcopy(nc))
                nc_neighbour[end].corner[dir]+=nc_neighbour[end].size[dir]
            end
            if issingchange(fcvals[indneg[dir]],Nf,Nc)
                push!(nc_neighbour,deepcopy(nc))
                nc_neighbour[end].corner[dir]-=nc_neighbour[end].size[dir]
            end
        end
    end
    #-------direcational neightbours----------------
    filter!(nc->!(any(nc.corner .< 1 ) || any((nc.corner+nc.size) .> [length.(mdbm.axes)...])),nc_neighbour)#remove the overhanging ncubes
    sort!(nc_neighbour; alg=QuickSort)
    # unique!(nc_neighbour)#TODO: ez miért nem jó, pedig definiálva van az "isequal"
    #unique!(a->[a.corner,a.size], nc_neighbour) #TODO: ez csak a Julia1.1 felett van!!!
    nc_neighbour=unique(a->[a.corner,a.size], nc_neighbour)
    return nc_neighbour
end

function checkneighbour!(mdbm::MDBM_Problem{N,Nf,Nc};interpolationorder::Int=0,maxiteration::Int=0) where N where Nf where Nc#only for unite size cubes
#mTODO: mivan azzal a kockával, aki szomszédként leellenőriztünk, de nem tartalmazott,... majd a követés során újra mellé kerülünk -> így az kétszer lesz ellőnőrizve (ez a Matlabban is rosssz)

if isempty(mdbm.ncubes)
    println("There is no bracketing n-cubes to check!")
else
    newbracketinncubes=mdbm.ncubes
    numberofiteration=0;
    while !isempty(newbracketinncubes) && (maxiteration==0 ? true : numberofiteration<maxiteration)
        numberofiteration+=1
        ncubes2check=generateneighbours(newbracketinncubes,mdbm)

        deleteat!(ncubes2check,is_sorted_in_sorted(ncubes2check,mdbm.ncubes))#delete the ones which is already presented

        _interpolate!(ncubes2check,mdbm, Val{interpolationorder})#remove the non-bracketing, only proper new bracketing cubes remained

        newbracketinncubes=deepcopy(ncubes2check)#TODO: kell a deepcopy?
        append!(mdbm.ncubes,deepcopy(ncubes2check))#TODO: kell a deepcopy?
        sort!(mdbm.ncubes; alg=QuickSort)
    end
end
end

function index_sorted_in_sorted(a::AbstractVector, b::AbstractVector)::Array{Int64,1}
    # index of a[i] in b (0 if not present)
    # a and b must be sorted
    containingindex=zeros(Int64,length(a))
    startindex=1;
    for ind2check in 1:length(a)
        detectedrange=searchsorted(b[startindex:end], a[ind2check])
        startindex=max(detectedrange.stop,detectedrange.start)+startindex-1
        containingindex[ind2check]=isempty(detectedrange) ? 0 : startindex
        if startindex>length(b)
            break
        end
    end
    return containingindex
end

function is_sorted_in_sorted(a::AbstractVector, b::AbstractVector)::Array{Bool,1}
    # is a[i] in b
    # a and b must be sorted
    iscontained=falses(length(a))
    startindex=1;
    for ind2check in 1:length(a)
        detectedrange=searchsorted(b[startindex:end], a[ind2check])
        iscontained[ind2check]=!isempty(detectedrange)
        startindex=max(detectedrange.stop,detectedrange.start)+startindex-1
        if startindex>length(b)
            break
        end
    end
    return iscontained
end




function DTconnect(mdbm::MDBM_Problem{N,Nf,Nc}) where N where Nf where Nc
    #---------- line connection (no corener neighbour is needed) --------------
    DT1=Array{Tuple{Int64,Int64}}(undef,0)
    for inc in 1:length(mdbm.ncubes)
        ncneigh=generateneighbours([mdbm.ncubes[inc]],mdbm)
        indinresults=index_sorted_in_sorted(ncneigh,mdbm.ncubes)#delete the ones which is already presented
        append!(DT1,[(inc,x) for x in indinresults if x!=0])
        # _interpolate!(ncubes2check,mdbm, Val{interpolationorder})#remove the non-bracketing, only proper new bracketing cubes remained
        #
        # newbracketinncubes=deepcopy(ncubes2check)#TODO: kell a deepcopy?
        # append!(mdbm.ncubes,deepcopy(ncubes2check))#TODO: kell a deepcopy?
    end
    filter!(d->(d[2]>d[1]),DT1)#TODO: eleve csak ez egyik irányban levő szomszédokat kellene keresni!!!
    sort!(DT1; alg=QuickSort)
    unique!(DT1)
    return DT1
end


function triangulation(DT1::Array{Tuple{Int64,Int64}})::Array{Tuple{Int64,Int64,Int64}}

    #DT=sort!()[DT1;[(d[2],d[1]) for d in DT1]]);
    DT=sort!([DT1;[d[[2,1]] for d in DT1]])#both direction of line connection is necessay!

    #L=[filter(d->d[1]==i,DT) for i in 1:length(mymdbm.ncubes)]
    L=[
    [dd[2] for dd in filter(d->d[1]==i,DT)] for i in 1:length(mymdbm.ncubes)]

    DT4=Array{Tuple{Int64,Int64,Int64,Int64}}(undef,0)#quadratic patch
    for i in 1:size(L,1)
        for j in L[i]#filter((L{i}>i) for i in )
            if i>j#i must be the larges value to remove the repetition of the some surface segment
                for k in L[j]#(L{j}>i) it seems to be slower
                    if i>k# there is no backstep, and i must be the largest   (Equivalent: %((k~=i) && (i>k))%back step is not allowed
                        for m in L[k]
                            if ((m!=j) && (j>m)) #&& (i>m) #back step is not allowed, and i must be the largest
                                if any(i.==L[m])
                                    push!(DT4,(i,j,k,m))
                                end
                            end
                        end
                    end
                end
            end
        end
    end

    return [[dt[[1,2,3]] for dt in DT4];[dt[[3,1,4]] for dt in DT4]]#DT2 triangular patch from the quadratic patch

end
