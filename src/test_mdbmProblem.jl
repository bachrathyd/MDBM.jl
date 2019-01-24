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

# ax1=Axis(1:5.0,"x")
# ax2=Axis(10:10:80.0,"y")
# ax3=Axis(1:10.0,"z")

mdbmaxes=[ax1,ax2,ax3]

function f(x,y,z)
    [x*x+y*y+z*z-2.0*2.0,x+y+z*2.0]
end

function c(x,y,z)
    [y,x+1.5]
end

# [mdbm.f,mdbm.c].(1.2,2.3,4.5)
# mdbm.c.memoryacc[1]

axdoubling!.(mdbmaxes)
# println(prod([length(mdbmaxes[i].ticks) for i in 1:3]))
@time mdbm=MDBM_Problem(f,mdbmaxes)
@time mdbm=MDBM_Problem(f,mdbmaxes,constraint=c)
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
x,y,z=getinterpolatedpoint(mdbm)
scatter(x,y,size = (2700, 2700))



# szemét--------------------------------

println("<<<<<<<<<<<<<<<<<<")
nc=mdbm.ncubes[2]
@code_warntype getcornerval2(nc,mdbm)
# mdbm.ncubes[2].size[:]=mdbm.ncubes[2].size.÷2
println(".................----------------------<<<<<<<<<<<<<<<<<<")
@code_warntype getcornerval(mdbm)
# @code_llvm getcornerval(mdbm)


mdbm.axes
println(".................----------------------<<<<<<<<<<<<<<<<<<")
@code_warntype mdbm.c(1.2,3.2,7.2)

@code_warntype mdbm.axes(11,11,1)

cf=mdbm.axes(11,1,1)
@code_warntype mdbm.f(cf...)
@code_warntype allcorner(nc)

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









# ax1=Axis([[0.0,1.0],[1.0,1.1]],"x")
# ax2=Axis(-7:1:7.0,"y")
#
# mdbmaxes=[ax1,ax2]
#
# function f(x,y)
#     [x*[1,y]]
# end
# interporder=1;#0
#
# @time mdbm=MDBM_Problem(f,mdbmaxes)
# @time interpolate!(mdbm,interpolationorder=interporder)
# for kstep=1:8
#     @time refine!(mdbm)
#     @time interpolate!(mdbm,interpolationorder=interporder)
# end
# x,y=getinterpolatedpoint(mdbm)
# println(length(x))
# scatter(x,y,size = (2700, 2700))
