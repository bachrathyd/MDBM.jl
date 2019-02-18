

function axdoubling!(ax::Axis)
        sort!(append!(ax.ticks, ax.ticks[1:end - 1] + diff(ax.ticks) / 2); alg=QuickSort)
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

function corner(nc::NCube{IT,FT,N},T01)::Vector{MVector} where IT where FT where N
    [nc.corner .+ nc.size .* T for T in T01]
end


function corner(ncubes::Vector{NCube{IT,FT,N}},T01)::Vector{Vector{MVector}} where IT where FT where N
    [corner(nc,T01) for nc in ncubes]
    # [nc.corner .+ nc.size .* T for nc in ncubes for T in T01]
end

function corner(mdbm::MDBM_Problem{N,Nf,Nc})::Vector{Vector{SVector}} where N where Nf where Nc
    [corner(nc,mdbm.T01) for nc in mdbm.ncubes]
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


function T01maker(valk::Val{kdim}) where {kdim}
        SVector{2^kdim}([
        SVector{kdim}([isodd(x÷(2^y)) for y in 0:(kdim-1)])
         for  x in 0:(2^kdim-1)])
end

function T11pinvmaker(valk::Val{kdim}) where {kdim}
    T11pinv=([isodd(x÷(2^y)) for y in 0:kdim , x in 0:(2^kdim-1)]*2.0 .-1.0)/(2^kdim)
    SMatrix{kdim+1,2^kdim}(T11pinv)
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

function _interpolate!(ncubes::Vector{NCube{IT,FT,N}},mdbm::MDBM_Problem{N,Nf,Nc},::Type{Val{0}}) where IT where FT where N where Nf where Nc
    isbracketing=map(nc->issingchange(getcornerval(nc,mdbm),Nf,Nc),ncubes) #TODO: parallelize
    deleteat!(ncubes,.!isbracketing)#removeing the non-bracketing ncubes

    for nc in ncubes
        nc.posinterp[:].=zero(FT)
    end
    return nothing
end

function _interpolate!(ncubes::Vector{NCube{IT,FT,N}},mdbm::MDBM_Problem{N,Nf,Nc},::Type{Val{1}}) where IT where FT where N where Nf where Nc

    for nc in ncubes
        FunTupleVector=getcornerval(nc,mdbm)

        if all([
            any((c)->!isless(c[2][fi],zero(c[2][fi])),FunTupleVector)
            for fi in 1:length(FunTupleVector[1][2])
            ])# check the constraint: do wh have to compute at all?!?

            As = zero(MVector{Nf+Nc,FT})
            ns = zero(MMatrix{N,Nf+Nc,FT})

            #for f---------------------
            for kf=1:Nf#length(FunTupleVector[1][1])
                solloc=mdbm.T11pinv*[FunTupleVector[kcubecorner][1][kf] for kcubecorner=1:length(FunTupleVector)]
                As[kf]=solloc[end];
                ns[:,kf].=solloc[1:end-1];
            end

            #for c---if needed: that is, it is close to the boundary--------
            activeCostraint=0;
            for kf=1:Nc#length(FunTupleVector[1][2])
                #TODO: mivan a többszörös C teljesülése esetén!?!??!# ISSUE
                if any((c)->!isless(zero(c[2][kf]),c[2][kf]),FunTupleVector) && length(As)<N  # use a constraint till it reduces the dimension to zero (point) and no further
                    solloc=mdbm.T11pinv*[FunTupleVector[kcubecorner][2][kf] for kcubecorner=1:length(FunTupleVector)]
                    activeCostraint+=1;
                    As[Nf+activeCostraint]=solloc[end];
                    ns[Nf+activeCostraint,kf].=solloc[1:end-1];
                end
            end

            if (Nf+activeCostraint)==0
                nc.posinterp[:] .=zero(FT)
            else
                nc.posinterp[:] .= transpose(ns[:,1:(Nf+activeCostraint)]) \As[1:(Nf+activeCostraint)];#TEST
            end
            #for c---------------------
        else
            nc.posinterp[:] .=1000.0;#put it outside the cube!
        end
    end

    #TODO: what if it falls outside of the n-cube, it should be removed ->what shall I do with the bracketing cubes?
    # Let the user define it
     # filter!((nc)->sum((abs.(nc.posinterp)).^10.0)<(1.5 ^10.0),mdbm.ncubes)#1e6~=(2.0 ^20.0)
     normp=10.0
     ncubetolerance=0.5
     filter!((nc)->norm(nc.posinterp,normp)<1.0+ncubetolerance,mdbm.ncubes)
     #filter!((nc)->!any(isnan.(nc.posinterp)),mdbm.ncubes)

    return nothing
end

function _interpolate!(ncubes::Vector{NCube{IT,FT,N}},mdbm::MDBM_Problem{N,Nf,Nc},::Type{Val{Ninterp}}) where IT where FT where N where Ninterp where Nf where Nc
    error("order $(Ninterp) interpolation is not supperted (yet)")
end

function interpolate!(mdbm::MDBM_Problem{N,Nf,Nc};interpolationorder::Int=1) where N where Nf where Nc
    _interpolate!(mdbm.ncubes,mdbm, Val{interpolationorder})
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
    for dir in directions
        for nc in ncubes
            nc.size[dir]/=2
        end
        NumofNCubes=length(ncubes)
        append!(ncubes,deepcopy(ncubes))
        for nc in ncubes[1:NumofNCubes]
            nc.corner[dir]=nc.corner[dir]+nc.size[dir]
        end
    end
    sort!(ncubes; alg=QuickSort)
    return nothing
end




function getinterpolatedsolution(mdbm::MDBM_Problem{N,Nf,Nc}) where N where Nf where Nc
[
    [
        (typeof(mdbm.axes[i].ticks).parameters[1])(
        (mdbm.axes[i].ticks[nc.corner[i]]*(1.0-(nc.posinterp[i]+1.0)/2.0)+
        mdbm.axes[i].ticks[nc.corner[i]+1]*((nc.posinterp[i]+1.0)/2.0))
        )
    for nc in mdbm.ncubes]#::Vector{typeof(mdbm.axes[i].ticks).parameters[1]}
for i in 1:length(mdbm.axes)]
end
function getinterpolatedsolution(nc::NCube{IT,FT,N},mdbm::MDBM_Problem{N,Nf,Nc}) where IT where FT where N where Nf where Nc
[
        (typeof(mdbm.axes[i].ticks).parameters[1])(
        (mdbm.axes[i].ticks[nc.corner[i]]*(1.0-(nc.posinterp[i]+1.0)/2.0)+
        mdbm.axes[i].ticks[nc.corner[i]+1]*((nc.posinterp[i]+1.0)/2.0))
        )#::Vector{typeof(mdbm.axes[i].ticks).parameters[1]}
for i in 1:length(mdbm.axes)]
end

function getevaluatedpoints(mdbm::MDBM_Problem{N,Nf,Nc}) where N where Nf where Nc
    [[x.callargs[i] for x in mdbm.fc.fvalarg] for i in 1:N]
end
function getevaluatedfunctionvalues(mdbm::MDBM_Problem{N,Nf,Nc}) where N where Nf where Nc
    [x.funval for x in mdbm.fc.fvalarg]
end
function getevaluatedconstraintvalues(mdbm::MDBM_Problem{N,Nf,Nc}) where N where Nf where Nc
    [x.cval for x in mdbm.fc.fvalarg]
end

# function generateneighbours(ncubes::Vector{NCube},mdbm::MDBM_Problem{N,Nf,Nc}) where N  where Nf where Nc
function generateneighbours(ncubes::Vector{NCube{IT,FT,N}},mdbm::MDBM_Problem{N,Nf,Nc}) where IT where FT where N  where Nf where Nc
    #TODO: Let the user select the method
    #-------all faceing/cornering neightbours----------------
    # #neighbourind=1:2^length(mdbm.axes) #cornering neightbours also - unnecessary
    # neighbourind=2 .^(0:(length(mdbm.axes)-1)) .+ 1 #neighbour only on the side
    # T101=[-mdbm.T01[neighbourind]...,mdbm.T01[neighbourind]...]
    #
    # nc_neighbour = Array{typeof(mdbm.ncubes[1])}(undef,0)
    # NumofNCubes=length(ncubes)
    # for iT in 1:length(T101)
    #     append!(nc_neighbour,deepcopy(ncubes))
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
    unique!(nc_neighbour)
    return nc_neighbour
end

function checkneighbour!(mdbm::MDBM_Problem{N,Nf,Nc};interpolationorder::Int=0,maxiteration::Int=0) where N where Nf where Nc#only for unite size cubes

if isempty(mdbm.ncubes)
    println("There is no bracketing n-cubes to check!")
else
    ncubes2check=mdbm.ncubes
    numberofiteration=0;
    while !isempty(ncubes2check) && (maxiteration==0 ? true : numberofiteration<maxiteration)
        numberofiteration+=1
        ncubes2check=generateneighbours(ncubes2check,mdbm)

        deleteat!(ncubes2check,is_sorted_in_sorted(ncubes2check,mdbm.ncubes))#delete the ones which is already presented
        _interpolate!(ncubes2check,mdbm, Val{interpolationorder})#remove the non-bracketing, only proper new bracketing cubes remained

        append!(mdbm.ncubes,deepcopy(ncubes2check))
        sort!(mdbm.ncubes; alg=QuickSort)
    end
end
end



function connect(mdbm::MDBM_Problem{N,Nf,Nc}) where N where Nf where Nc
    #---------- line connection (no corener neighbour is needed) --------------
    DT1=Array{Tuple{Int64,Int64}}(undef,0)
    for inc in 1:length(mdbm.ncubes)
        ncneigh=generateneighbours([mdbm.ncubes[inc]],mdbm)
        indinresults=index_sorted_in_sorted(ncneigh,mdbm.ncubes)#delete the ones which is already presented
        append!(DT1,[(inc,x) for x in indinresults if x!=0])
    end
    filter!(d->(d[2]>d[1]),DT1)#TODO: eleve csak ez egyik irányban levő szomszédokat kellene keresni!!!
    sort!(DT1; alg=QuickSort)
    unique!(DT1)
    return DT1
end


function triangulation(DT1::Array{Tuple{Int64,Int64}})::Array{Tuple{Int64,Int64,Int64}}

    #DT=sort!()[DT1;[(d[2],d[1]) for d in DT1]]);
    DT=sort!([DT1;[d[[2,1]] for d in DT1]])#both direction of line connection is necessay!

    #L=[filter(d->d[1]==i,DT) for i in 1:length(mdbm.ncubes)]
    L=[
    [dd[2] for dd in filter(d->d[1]==i,DT)] for i in 1:maximum(DT1)[2]]

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



function solve!(mdbm::MDBM_Problem{N,Nf,Nc},iteration::Int;interpolationorder::Int=1) where N where Nf where Nc
    interpolate!(mdbm,interpolationorder=interpolationorder)
    for k=1:iteration
        refine!(mdbm)
        interpolate!(mdbm,interpolationorder=interpolationorder)
    end
    checkneighbour!(mdbm)
    interpolate!(mdbm,interpolationorder=interpolationorder)
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
