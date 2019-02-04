include("MDBM__types.jl")
# ~~~~~~~~~~~~~~~~~~~~~~~
ax1=Axis(1:20)
ax2=Axis(1:20)
ax3=Axis(1:-1:-20)

mdbmaxes=[ax1,ax2,ax3]

foo1=(mdbmaxes)->[NCube([x[1]...],ones(Int64,length(x[1])),zeros(Float64,length(x[1])))
    for x in Iterators.partition(Iterators.product((x->eachindex(x.ticks)).(mdbmaxes)...),1)]
foo2=(mdbmaxes)->[NCube(Int16.([x...]),ones(Int16,length(x)),zeros(Float64,length(x)))
    for x in Iterators.product((x->eachindex(x.ticks)).(mdbmaxes)...)][:]
foo3=(mdbmaxes)->[NCube(Int16.([x...]),ones(Int16,length(x)),zeros(Float64,length(x)))
                 for x in Iterators.product((x->eachindex(x.ticks)).(mdbmaxes)...)]|>vec
println("-<<<<<<<<>>>>>>>>>>>--")
@time foo1(mdbmaxes)
@time foo2(mdbmaxes)
@time foo3(mdbmaxes)
println("-................")

# ~~~~~~~~~~~~~~~~~~~~~~~


include("MDBM__types.jl")

ax1=Axis(-5:5.0,"x")
ax2=Axis(Float16.(-7:1:7.0),"y")
ax3=Axis(-3:3.0,"z")

mdbmaxes=[ax1,ax2,ax3]
function f(x,y,z)
    [x*x+y*y+z*z-2.0*2.0,x+y+z*2.0]
end
function c(x,y,z)
    [y,x+1.5]
end

@time mdbm=MDBM_Problem(f,mdbmaxes,constraint=c)



ax1=Axis(-5:5.0,"x")
ax2=Axis(Float16.(-7:1:7.0),"y")

mdbmaxes=[ax1,ax2]
function f(x,y)
    [x-y,y-1]
end

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

for k=1:4
println("-................")
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


    [
        mdbm.axes[i].ticks[nc.corner[i]]*(1-(nc.posinterp[i]+1)/2)+
        mdbm.axes[i].ticks[nc.corner[i]+1]*((nc.posinterp[i]+1)/2)

for i in 1:length(mdbm.axes)]

# szemét--------------------------------

println("<<<<<<<<<<<<<<<<<<")
nc=mdbm.ncubes[1]



@code_warntype getcornerval(nc,mdbm)
# mdbm.ncubes[2].size[:]=mdbm.ncubes[2].size.÷2
println(".................----------------------<<<<<<<<<<<<<<<<<<")
@code_warntype getcornerval2(nc,mdbm)


# @code_llvm getcornerval(mdbm)

nc.corner::Vector{IT}.+((mdbm.T01).*(nc.size::Vector{IT}))
mdbm.axes
println(".................----------------------<<<<<<<<<<<<<<<<<<")
@code_warntype mdbm.c(1.2,3.2,7.2)

@code_warntype mdbm.axes(11,11,1)
@code_warntype mdbm.axes([11,11,1]...)
T=allcorner(nc,mdbm.T01)

Ndim=3
T01=[[isodd(x÷(2^y)) for y in 0:(Ndim-1)] for x in 0:(2^Ndim-1) ]

Base.:*.(T01,Ref(nc.size))

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
