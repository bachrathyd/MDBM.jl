#Light graph test


using LightGraphs


# myDT1=connect(mymdbm)

g = SimpleGraph(length(mymdbm.ncubes))
dg = SimpleDiGraph(length(mymdbm.ncubes));

for dt in myDT1
add_edge!(g, dt[1], dt[2]);
add_edge!(dg, dt[2], dt[1]);
add_edge!(dg, dt[1], dt[2]);
end




fig = figure(13);clf()

a=connected_components(g)

for i in 1:length(a)

    Ps=getinterpolatedsolution(mymdbm.ncubes[a[i]],mymdbm)
    # P2=getinterpolatedsolution(mymdbm.ncubes[dt[2]],mymdbm)
    plot3D(Ps[1],Ps[2],Ps[3],marker=".")
end

is_cyclic(g)

# a=simplecycleslength(g, 10^6)
# simplecycles_hawick_james(dg)
