using StaticArrays
using LinearAlgebra


#
#println("Original faces (3-cubes):")
#for (i, face) in enumerate(faces)
#    println("Face $(i), fixed dim: ", face_fixed_dims[i])
#end
#
#println("\nSub-faces of first face (2-cubes):")
#sub_faces, all_fixed_dims = generate_sub_faces(faces[1], all_fixed_dims[1])
#for (i, sub_face) in enumerate(sub_faces)
#    println("Sub-face $(i), all fixed dims: ", all_fixed_dims[i])
#end
#
#println("\nEdges (1-cubes) from first sub-face:")
#edges, all_edge_fixed_dims, = generate_sub_faces(sub_faces[1], all_fixed_dims[1])
#for (i, edge) in enumerate(edges)
#    println("Edge $(i), all fixed dims: ", all_edge_fixed_dims[i])
#end
#




#fighanle = scatter()
fighanle = scatter!()

N = 4
T = Float64
pathall = Vector{SVector{N,T}}()
for nc in mymdbm.ncubes
#nc= mymdbm.ncubes[45]
    nccorner = nc.corner
    ncsize = nc.size

    fixed_dims = Vector{Tuple{Int,Bool}}(undef, 0)
    nc_template = mymdbm.T01


    posall_tree = MDBM.PositionTree(NaN * Vector{Float64}(undef, N))
    MDBM.interpsubcubesolution!(posall_tree, [nc_template], [fixed_dims], nccorner, ncsize, mymdbm)
    #    
#    @show posall_tree
    path2points = MDBM.extract_paths(posall_tree)

    #[plot!(collect(eachrow(reduce(hcat, pathloc)))...) for pathloc in path2points]
    #plot!() 
    
    #pathloc=path2points[1]
    for pathloc in path2points
        if length(pathloc)>2
        #append!(pathall, pathloc)
        append!(pathall, vcat(pathloc, [pathloc[2]]))
        end
    end

end


plot(collect(eachrow(reduce(hcat, pathall)))[1:3]..., lw=5)
plot!(collect(eachrow(reduce(hcat, pathall)))[1:2]..., lw=5)
#scatter!(collect(eachrow(reduce(hcat, pathall)))..., ms=15, label="interpolated solution")