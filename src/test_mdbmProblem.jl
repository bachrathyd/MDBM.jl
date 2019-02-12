# include("MDBM__types.jl")
# # ~~~~~~~~~~~~~~~~~~~~~~~
# ax1=Axis(1:20)
# ax2=Axis(1:20)
# ax3=Axis(1:-1:-20)
#
# mdbmaxes=[ax1,ax2,ax3]
#
# foo1=(mdbmaxes)->[NCube([x[1]...],ones(Int64,length(x[1])),zeros(Float64,length(x[1])))
#     for x in Iterators.partition(Iterators.product((x->eachindex(x.ticks)).(mdbmaxes)...),1)]
# foo2=(mdbmaxes)->[NCube(Int16.([x...]),ones(Int16,length(x)),zeros(Float64,length(x)))
#     for x in Iterators.product((x->eachindex(x.ticks)).(mdbmaxes)...)][:]
# foo3=(mdbmaxes)->[NCube(Int16.([x...]),ones(Int16,length(x)),zeros(Float64,length(x)))
#                  for x in Iterators.product((x->eachindex(x.ticks)).(mdbmaxes)...)]|>vec
# println("-<<<<<<<<>>>>>>>>>>>--")
# @time foo1(mdbmaxes)
# @time foo2(mdbmaxes)
# @time foo3(mdbmaxes)
# println("-................")
#
# # ~~~~~~~~~~~~~~~~~~~~~~~


include("MDBM__types.jl")

println("-................")
#TDOO: WTF
# ax1=Axis((-5,5),"x")
# ax1=Axis("asdf","ppp")
ax1=Axis([-5,0,3,5],"x")
ax2=Axis(Float16,(-7:1:7.0))
ax3=Axis(-3:3.0,"z")
#TODO: példában megmutatni, de szólni, hogy nem lehet int (interpolálni kell tudni)
# ax3=Axis([[1,5],[2,3.0]],"z")
# ax3=Axis(Int64,-3:4,"z")

mdbmaxes=[ax1,ax2,ax3]

function f(x,y,z)
    #[x*x+y*y+z*z-2.0*2.0]

    pN=10.85
    [abs(x)^pN+abs(y)^pN+abs(z)^pN-2.0^pN]

    #[x*x+y*y+z*z-2.0*2.0,(x+y+z*2.0)]#,y<0,(x+1.5)<0]
    #[x*x+y*y+z*z-2.0*2.0,(x+y+z*2.0),y<0,(x+1.5)<0] #így is működik
end
function c(x,y,z)
    minimum([y,x+1.5,z])
end

# @time mymdbm=MDBM_Problem(f,mdbmaxes)
# @time mymdbm=MDBM_Problem(f,mdbmaxes,constraint=c,memoization=false)
@time mymdbm=MDBM_Problem(f,mdbmaxes,constraint=c)

# ax1=Axis(-5:5.0,"x")
# ax2=Axis(Float16.(-7:1:7.0),"y")
#
# mdbmaxes=[ax1,ax2]
# function f(x,y)
#     [x-y,y-1]
# end
#
# @time mdbm=MDBM_Problem(f,mdbmaxes)


# [mdbm.f,mdbm.c].(1.2,2.3,4.5)
# mdbm.c.memoryacc[1]

# println(prod([length(mdbmaxes[i].ticks) for i in 1:3]))
# @time mdbm=MDBM_Problem(f,mdbmaxes)
# @time mdbm=MDBM_Problem(f,mdbmaxes,constraint=c)
# @time mdbm=MDBM_Problem(f,mdbmaxes,constraint=c,memoization=false)
# @time mdbm=MDBM_Problem(f,mdbmaxes)
# @time interpolate!(mdbm,interpolationorder=0)
InterporderN=1
@time interpolate!(mymdbm,interpolationorder=InterporderN)


@time for k=1:5
refine!(mymdbm)
interpolate!(mymdbm,interpolationorder=InterporderN)
println("-...==========......")
end

using Plots
gr()

x,y,z=getinterpolatedpoint(mymdbm)
println("pontins: ", length(mymdbm.ncubes))
println("function evaluation and memo: ",[length(mymdbm.fc.fvalarg),mymdbm.fc.memoryacc[1]])
scatter(x,y,z)
println(length(mymdbm.fc.fvalarg)*5)
println("datapoints: ", length(mymdbm.ncubes)*9+length(mymdbm.fc.fvalarg)*5)

#--------------- szemét--------------------------------
# @time mdbm.fc(1.0,2.0,2.0)
# @time mdbm.fc(1.0,2.0,rand(Float64,1)...)
# @time mdbm.fc(1.0,2.0,3.0)]#

aa=mymdbm.ncubes[2]
bb=deepcopy(mymdbm.ncubes[2])
bb.posinterp[:]+=10.00;
aa>=bb
isequal(aa,bb)
aa==bb
@code_warntype  compa(aa,bb)
all([a.corner==b.corner,a.size==b.size])# @code_warntype mdbm.fc(1.0,2.0,3.0)
# @code_warntype mdbm.fc(1.0,2.0,rand(Float64,1)...)
filter!((nc)->any(isnan.(nc.posinterp)),mymdbm.ncubes)
nc
println("<<<<<<<<<<<<<<<<<<")
nc=mymdbm.ncubes[1]

mdbm=mymdbm
T101=mdbm.T01
newbracketinncubes =mdbm.ncubes
ncubes2check = Array{typeof(mdbm.ncubes[1])}(undef,0)
NumofNCubes=length(newbracketinncubes)
for iT in 1:length(T101)
    append!(ncubes2check,deepcopy(newbracketinncubes))#TODO: ez így kell csinálni?
    for nci in ((1+NumofNCubes*(iT-1)):(NumofNCubes+NumofNCubes*(iT-1)))
        ncubes2check[nci].corner[:]=ncubes2check[nci].corner+T101[iT].*ncubes2check[nci].size
            # ncubes2check[nci].corner[:]=ncubes2check[nci].corner+T.*ncubes2check[nci].size
    end
end


        sort!(ncubes2check; alg=QuickSort)
        # unique!(ncubes2check)#TODO: ez miért nem jó, pedig definiálva van az "isequal"
        #unique!(a->[a.corner,a.size], ncubes2check) #TODO: ez csak a Julia1.1 felett van!!!
        ncubes2check=unique(a->[a.corner,a.size], ncubes2check)
        is_sorted_in_sorted(ncubes2check,mdbm.ncubes)
        sum(is_sorted_in_sorted(ncubes2check,mdbm.ncubes))

deleteat!(ncubes2check,is_sorted_in_sorted(ncubes2check,mdbm.ncubes))

_interpolate!(ncubes2check,mdbm, Val{1})

    findin(ncubes2check,mdbm.ncubes)
newneighbour = cubinds[indexin_sorted_iszero(mdbm.ncubes, ncubes2check)]
findin([1,43,3,6,7],collect(2:14))
findin(ncubes2check,mdbm.ncubes)
searchsortedfirst()

b = [1,2,5,7,9]
a=[2,3,4]

nc2chk=ncubes2check[3000]

println("-----%%%%%%%%%%%%%%%%%-------")
is_sorted_in_sorted([1,2,3,4,5,7,7,7,7,10,12,13,15],[1,2,5,7,7,7,7,9])


println("-zzzzzzzzzzzzzzzzzzzzzzzzzzzzz----%%%%%%%%%%%%%%%%%-------")
is_sorted_in_sorted([1,1,2,3,5,7,7,7,7,10,12,13,15],[1,2,5,7,7,7,7,9])


Base.Sort.searchsorted(a,b[1])
mdbm.ncubes[2]>=ncubes2check[2]
filter(x->,ncubes2check)

aaa=[any(nc->nc2chk==nc, mdbm.ncubes) for nc2chk in ncubes2check]
any(aaa)
mdbm=mymdb
ű-., m
interpolationorder=1

allcorner(nc,[-mymdbm.T01[2:end],mymdbm.T01[2.^(0:2)+1]])

doubling!(mdbm,[1,1,2])
refinecubes!(mdbm.ncubes,[1])

refine!(mdbm)

typeof(allcorner(nc,mdbm.T01))
@code_warntype allcorner(nc,mdbm.T01)
@code_warntype allcorner(mdbm.ncubes,mdbm.T01)

@code_warntype getcornerval(mdbm.ncubes[2],mdbm)
@code_warntype getcornerval(mdbm)

@code_warntype map((mdbm.fc) ∘ (mdbm.axes),allcorner(nc,mdbm.T01))

@code_warntype allcorner(mdbm)
@code_warntype ((allcorner(mdbm.ncubes,mdbm.T01)))
T=T101[1]
@code_warntype (mdbm.axes.(allcorner(nc,mdbm.T01)))
@code_warntype map((mdbm.axes),Base.Iterators.flatten(allcorner(mdbm)))
typeof(map((mdbm.axes),Base.Iterators.flatten(allcorner(mdbm))))
(mdbm.fc).(map((mdbm.axes),Base.Iterators.flatten(allcorner(mdbm))))
map((mdbm.fc) ∘ (mdbm.axes),Base.Iterators.flatten(allcorner(mdbm)))
@code_warntype (mdbm.fc).(map((mdbm.axes),Base.Iterators.flatten(allcorner(mdbm))))
@code_warntype map((mdbm.fc) ∘ (mdbm.axes),Base.Iterators.flatten(allcorner(mdbm)))
@time (mdbm.fc).(map((mdbm.axes),Base.Iterators.flatten(allcorner(mdbm))))
@time map((mdbm.fc) ∘ (mdbm.axes),Base.Iterators.flatten(allcorner(mdbm)))

(mdbm.fc). ((mdbm.axes).(allcorner(mdbm.ncubes[1],mdbm.T01)))

@time

mdbm.fc(111.7,Float16(1.2),31.3)
[mdbm.fc.memoryacc[1],length(mdbm.fc.fvalarg)]

mdbm.fc
mdbm.fc((2.0,Float16(1.6),3.3))
mdbm.fc([(2.0,Float16(1.6),3.3),(5.0,Float16(1.6),3.3),(2.0,Float16(1.6),3.3)])



location=searchsortedfirst(mdbm.fc.fvalarg,(1.7, 1.2, 3.3),lt=(x,y)->isless(x.callargs,y));



collect(Base.Iterators.flatten(allcorner(mdbm)))
@code_warntype (x->mdbm.fc(x...)).(map((mdbm.axes),Base.Iterators.flatten(allcorner(mdbm))))

cc=allcorner(nc,mdbm.T01)[1]
mdbm.axes.(cc)
typeof(mdbm.ncubes[1].corner)

mdbm.T01[1]
@code_warntype mdbm.axes(3,4,2)
println(".................----------------------<<<<<<<<<<<<<<<<<<")




map(mdbm.axes,[1,1,3])
@code_warntype getcornerval(mdbm)
@time _interpolate!(mdbm,Val{1})
typeof(mdbm.ncubes[1].corner[1])
# @code_llvm getcornerval(mdbm)

nc.corner::Vector{IT}.+((mdbm.T01).*(nc.size::Vector{IT}))
@code_warntype mdbm.axes(1,2,3)
@code_warntype mdbm.axes([1,2,3])
@code_warntype mdbm.axes((1,2,3))
println(".................----------------------<<<<<<<<<<<<<<<<<<")
@code_warntype mdbm.c(1.2,3.2,7.2)

@code_warntype mdbm.axes(11,11,1)
@code_warntype mdbm.axes([11,11,1]...)
T=allcorner(nc,mdbm.T01)
    fcvals=getcornerval(nc,mdbm)
    fcvals[1]
    kf=1

        solloc=mdbm.T11pinv*[fcvals[1][kcube][kf] for kcube=1:length(fcvals[1])]
        [fcvals[1][:][kf]]

Ndim=3
T01=[[isodd(x÷(2^y)) for y in 0:(Ndim-1)] for x in 0:(2^Ndim-1) ]

Base.:*.(T01,Ref(nc.size))
typeof(nc.corner)

Vector{Float32}(2)
T01[:].*nc.size

nc.corner::Vector{IT}.+(((nc.size::Vector{IT}*).()
# using Makie
# scene = Scene()
# scene = scatter(x, y, color = colors)


    # plot!(title = "New Title", xlabel = "New xlabel", ylabel = "New ylabel")
    # plot!(xlims = (0, 5.5), ylims = (-2.2, 6), xticks = 0:0.5:10, yticks = [0,1,5,10])
    #
    # # or using magic:
    # plot!(xaxis = ("mylabel", :log10, :flip))
    # xaxis!("mylabel", :log10, :flip)


#     @time getinterpolatedpoint(mdbm)
println("<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<")
@code_warntype getinterpolatedpoint(mdbm)

kdim=2

@code_warntype Tmaker2(Val(2))

typeof(Tmaker(2))

#zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz
axtest=Axis([Float32.([1.0,0.0]),Float32.([1.0,3.0])],"vec")
# axtest=Axis([0.0,3.0],"scal")
axdoubling!(axtest)
ax2=Axis(Float16.(-7:1:7.0),"y")

using LinearAlgebra

mdbmaxes=[axtest,ax2]
function f(x,y)
    (norm(x)-y)>0.0
end
# function f(x,y)
#     sqrt(1.0+x*x)-y
# end

@time mdbm=MDBM_Problem(f,mdbmaxes)
# [mdbm.f,mdbm.c].(1.2,2.3,4.5)
# mdbm.c.memoryacc[1]

# println(prod([length(mdbmaxes[i].ticks) for i in 1:3]))
# @time mdbm=MDBM_Problem(f,mdbmaxes)
# @time mdbm=MDBM_Problem(f,mdbmaxes,constraint=c)
# @time mdbm=MDBM_Problem(f,mdbmaxes,constraint=c,memoization=false)
# @time mdbm=MDBM_Problem(f,mdbmaxes)
# @time interpolate!(mdbm,interpolationorder=0)
@time interpolate!(mdbm,interpolationorder=1)

println("-................")
for k=1:4
@time refine!(mdbm)
@time interpolate!(mdbm,interpolationorder=1)
println("-...==========......")
end

using Plots
gr()
# plotlyjs()
# unicodeplots()
# x,y,z=getinterpolatedpoint(mdbm)
x,y=getinterpolatedpoint(mdbm)
scatter(x,y)

@code_warntype getinterpolatedpoint(mdbm)




# SH test
include("MDBM__types.jl")
using Plots
gr()

using LinearAlgebra
mymdbm=MDBM_Problem((x...)->norm([x...].+ 0.1,2.7)-2,[-2:2,-2:2,-2:2])

interpolate!(mymdbm,interpolationorder=1)
for k=1:2
    @time refine!(mymdbm)
    @time interpolate!(mymdbm,interpolationorder=1)
    println("------")
end
@time checkneighbour!(mymdbm)
