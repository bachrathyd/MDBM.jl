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

ax1=Axis(1:6.0,"x")
ax2=Axis(10:10:30.0,"y")
ax3=Axis([100,200.0],"z")

mdbmaxes=[ax1,ax2,ax3]

function f(x,y,z)
    [x,y]
end

function c(x,y,z)
    y
end


mdbm=MDBM_Problem(f,mdbmaxes)

mdbm=MDBM_Problem(f,a)


a=[1:6.0,10:10:30.0,[100,200.0]]
mdbm=MDBM_Problem((x,y,z)->Float64((x^2+y^2.0-1.0)),a)




ccc=[(mdbm.ncubes[end].corner[n],mdbm.ncubes[end].corner[n]+mdbm.ncubes[end].size[n]) for n in 1:length(mdbm.ncubes[end].corner)]

[(mdbm.f).((mdbm.axes).(x)...) for x in Iterators.product(ccc...)]


fvals=Vector
for nc in mdbm.ncubes
println(
 [
 (mdbm.f).((mdbm.axes).(x)...)
 for x in
 Iterators.product(
             [(nc.corner[n],nc.corner[n]+nc.size[n])
             for n in 1:length(nc.corner)]
     ...)
 ][:][:]
 )
end


map((x,y)->x.ticks[y], mdbm.axes,(1,2,1))
(mdbm.axes).(ccc)


(mdbm.ncubes[7].corner...)
(x->(x.corner).(mdbm.ncubes)
[(x->(x.corner...)).(mdbm.ncubes)]|>(mdbm.f)

mdbm.f(1.3,1.2)
