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
checkneighbour!(mymdbm)
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
nc=mymdbm.ncubes[1]
@code_warntype (corner(nc,mymdbm.T01))
