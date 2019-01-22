# MDBM__types.jl

struct MDBMcontainer{RT,AT}
    funval::RT
    callargs::AT
end

struct MemF{RT,AT}# <:Function
    f::Function
    fvalarg::Vector{MDBMcontainer{RT,AT}}#funval,callargs
end



#------------SH version - hidden 'Any'-------------
(memfun::MemF{RT,AT})(::Type{RT},args...) where {RT,AT} = memfun.f(args...)::RT

function (memfun::MemF{RT,AT})(args...) where {RT,AT}
    location=searchsortedfirst(memfun.fvalarg,args,lt=(x,y)->isless(x.callargs,y));
    if length(memfun.fvalarg)<location
        #println("ujraszamolas a végére")
        x=memfun(RT,args...);
        push!(memfun.fvalarg,MDBMcontainer{RT,AT}(x,args))
        return x
    elseif  memfun.fvalarg[location].callargs!=args
        #println("ujraszamolas közé illeszt")
        x=memfun(RT,args...);
        insert!(memfun.fvalarg, location, MDBMcontainer{RT,AT}(x,args))
        return x
    else
        #println("mar megvolt")
        return memfun.fvalarg[location].funval;
    end
end
#-------------------------

# function (memfun::MemF{RT,AT})(Vector{T}) where {T,RT,AT}
#     locations=indexin_sorted(memfun.fvalarg,args,lt=(x,y)->isless(x.callargs,y));
#     for x in locations if x!=0
#         x=memfun(RT,args...)::RT;
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
Base.setindex!(ax::Axis, X, inds...) = setindex!(ax.ticks, X, inds...)
Base.size(ax::Axis) = size(ax.ticks)


function Axis(a::Vector{T}) where T
    Axis(a, :unknown)
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



struct MDBM_Problem
    f::MemF
    c::MemF
    axes::Vector{Axis}
    ncubes::Vector{NCube}
    function MDBM_Problem(f,c,axes,ncubes::Vector{NCube{IT,FT}}) where IT where FT
        println(mdbmaxes)
        new(f,c,axes,[NCube(IT.([x...]),ones(IT,length(x)),zeros(FT,0)) for x in Iterators.product((x->1:(length(x.ticks)-1)).(mdbmaxes)...)][:])
    end
end
function (axes::Vector{Axis})(i::Vector)
    map((x,y)->x.ticks[y], mdbm.axes,i)
end

function (axes::Vector{Axis})(i...)
    map((x,y)->x.ticks[y], mdbm.axes,i)
end


# function MDBM_Problem{IT,FT}(f::Function, mdbmaxes::Vector{Axis};constraint::Function=(x...)->Float16(1.))
function MDBM_Problem(f::Function, mdbmaxes::Vector{<:Axis};constraint::Function=(x...)->Float16(1.))
    println("Most mdbmaxes::Vector{Axis} fut")
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


    fun=MemF(f,Array{MDBMcontainer{RTf,AT}}(undef, 0));
    cons=MemF(constraint,Array{MDBMcontainer{RTc,AT}}(undef, 0));


    MDBM_Problem(fun,cons,mdbmaxes,Vector{NCube{Int64,Float64}}(undef, 0))
end

# {AbstractArray{T,1} where T}
# function MDBM_Problem(f::Function, a::AbstractVector;constraint::Function=(x...)->Float16(1.))
function MDBM_Problem(f::Function, a::Vector{<:AbstractVector};constraint::Function=(x...)->Float16(1.))
    mdbmaxes=[Axis(ax) for ax in a]
    println("Most  a::Array{AbstractArray{T,1} where T,1} fut")
    argtypesofmyfunc=map(x->typeof(x).parameters[1], mdbmaxes);#Argument Type
    AT=Tuple{argtypesofmyfunc...};
    type_f=Base.return_types(f,AT)
    if length(type_f)==0
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

    fun=MemF(f,Array{MDBMcontainer{RTf,AT}}(undef, 0));
    cons=MemF(constraint,Array{MDBMcontainer{RTc,AT}}(undef, 0));

    MDBM_Problem(fun,cons,mdbmaxes,Vector{NCube{Int64,Float64}}(undef, 0))
end



function _interpolate!(mdbm,::Type{Val{0}})
    print("0 cuccc happening")
end

function _interpolate!(mdbm,::Type{Val{1}})
    print("1 cuccc happening")
    Ndim=length(mdbm.axes)

    axlocdimles = Array{Float64}(undef, 0)
    for kdim = 1:Ndim
        axlocdimles = hcat([axlocdimles;-ones(1, 2^(kdim - 1))], [axlocdimles;ones(1, 2^(kdim - 1))]);
    end
    TAn = hcat(-ones(2^3, 1), transpose(axlocdimles))
    #     TAn2=inv(TAn.'*TAn);
    #     TAtrafo=TAn2*TAn.';
    TAtrafo = (transpose(TAn) * TAn) \ transpose(TAn)
end

function _interpolate!(mdbm,::Type{Val{N}}) where N
    error("order $(N) interpolation is not supperted (yet)")
end

function interpolate!(mdbm::MDBM_Problem;interpolationorder::Int=1)
    _interpolate!(mdbm, Val{interpolationorder})
end

function getcornerval(mdbm::MDBM_Problem)

    for nc in mdbm.ncubes
     [
     (mdbm.f).((mdbm.axes).(x)...)
     for x in
     Iterators.product(
                 [(nc.corner[n],nc.corner[n]+nc.size[n])
                 for n in 1:length(nc.corner)]
         ...)
     ][:]
    end


end









#
#
#
#
# function interataion(mdbm::MDBM_Problem)(Nitertaion::Int=1)
#     for i in 1:Nitertaion
#
#     end
#
# end
