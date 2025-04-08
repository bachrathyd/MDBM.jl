using StaticArrays

"""
    generate_faces_indices(n)

Generates the faces of an n-dimensional cube.

# Arguments
- `n::Int`: Dimension of the cube.

# Returns
- `faces`: Vector containing the faces, each face represented by a vector of corners.
- `fixed_dims`: Vector containing tuples indicating the dimension and the side (0 or 1) fixed for each face.
"""
function generate_faces_indices(n)
    corners = [SVector{n}(digits(i, base=2, pad=n) .== 1) for i in 0:(2^n-1)]
    faces = Vector{Vector{SVector{n,Bool}}}(undef, 2 * n)
    fixed_dims = Vector{Tuple{Int,Bool}}(undef, 2 * n)

    idx = 1
    for dim in 1:n
        for side in [false, true]
            face = filter(c -> c[dim] == side, corners)
            faces[idx] = face
            fixed_dims[idx] = (dim, side)
            idx += 1
        end
    end

    return faces, fixed_dims
end

"""
    generate_sub_faces(face, fixed_dims)

Generates sub-faces of a given face by fixing additional dimensions.

# Arguments
- `face`: The current face represented by a vector of corners.
- `fixed_dims::Vector{Tuple{Int,Bool}}`: Dimensions already fixed in the current face along with their sides.

# Returns
- `sub_faces`: Vector containing sub-faces derived from the current face.
- `sub_fixed_dims`: Vector containing tuples indicating the new dimension and side fixed for each sub-face.
- `all_fixed_dims`: Vector containing tuples of all dimensions and sides fixed including previous steps.
"""
function generate_sub_faces(face, fixed_dims::Vector{Tuple{Int,Bool}})
    n = length(face[1])
    sub_faces = Vector{Vector{SVector{n,Bool}}}(undef, 0)
    sub_fixed_dims = Vector{Tuple{Int,Bool}}(undef, 0)
    fixed_dimensions_only = [fd[1] for fd in fixed_dims]
    free_dims = setdiff(1:n, fixed_dimensions_only)

    for dim in free_dims
        for side in [false, true]
            sub_face = filter(c -> c[dim] == side, face)
            push!(sub_faces, sub_face)
            push!(sub_fixed_dims, (dim, side))
        end
    end

    all_fixed_dims = [vcat(fixed_dims, [sf]) for sf in sub_fixed_dims]

    return sub_faces, all_fixed_dims#, sub_fixed_dims, 
end


function interpsubcubesolution(faces, fixed_dims, corner, axes, Nf, Nc, N)
    #(face, fixdim) = (faces[1], fixed_dims[1])
    posall_tree = Any[]

    for (face, fixdim) in zip(faces, fixed_dims)
        FunTupleVector = MDBM.getcornerval([nc.corner .+ nc.size .* T for T in face], mymdbm)

        fixed_dimensions_only = [fd[1] for fd in fixdim]
        free_dims = setdiff(1:N, fixed_dimensions_only)
        Nfree = length(free_dims)

        T11pinv = MDBM.T11pinvmaker(Val(length(free_dims)))



        posinterp = zeros(Float64, N)
        posinterp[getindex.(fixdim, 1)] .= -1 .+ 0.2 .* getindex.(fixdim, 2)
        posinterp[free_dims] .= MDBM.fit_hyperplane(FunTupleVector, Nfree, Nf, Nc, typeof(FunTupleVector[1][1][1]), T11pinv)
        posi = getinterpolatedsolution(posinterp, corner, axes)
        # scatter!([x0],[y0],[z0],ms=5)
        #x0,y0=getinterpolatedsolution(nc.posinterp, nc.corner, mymdbm.axes)
        #scatter!([x0],[y0],ms=5)

        #[plot!(-1 .+ 2 .* getindex.(f, 1), -1 .+ 2 .* getindex.(f, 2), -1 .+ 2 .* getindex.(f, 3)) for f in faces]
        # scatter!()

        normp = 50000.0
        ncubetolerance = 0.001
        # @show posinterp
        # @show norm(posinterp, normp) < 1.0 + ncubetolerance
        if norm(posinterp, normp) < 1.0 + ncubetolerance
            print("$Nfree ok:")
            @show face

            if Nfree > Nf
                edges, all_edge_fixed_dims, = generate_sub_faces(face, fixdim)
                posall_tree_sub = interpsubcubesolution(edges, all_edge_fixed_dims, corner, axes, Nf, Nc, N)
                push!(posall_tree, [posinterp, posall_tree_sub])
            end
            #x, y, z = getinterpolatedsolution(posinterp, nc.corner, mymdbm.axes)
            #plot!([x0, x], [y0, y], [z0, z], lw=4)
            #x,y=getinterpolatedsolution(posinterp, nc.corner, mymdbm.axes)
            #scatter!([x],[y],[z])
            #plot!([x0,x],[y0,y],lw=4)

        else
            println("$Nfree x")
            push!(posall_tree, [posinterp, []])
        end
    end
    return posall_tree
end


println("Original faces (3-cubes):")
for (i, face) in enumerate(faces)
    println("Face $(i), fixed dim: ", face_fixed_dims[i])
end

println("\nSub-faces of first face (2-cubes):")
sub_faces, all_fixed_dims = generate_sub_faces(faces[1], all_fixed_dims[1])
for (i, sub_face) in enumerate(sub_faces)
    println("Sub-face $(i), all fixed dims: ", all_fixed_dims[i])
end

println("\nEdges (1-cubes) from first sub-face:")
edges, all_edge_fixed_dims, = generate_sub_faces(sub_faces[1], all_fixed_dims[1])
for (i, edge) in enumerate(edges)
    println("Edge $(i), all fixed dims: ", all_edge_fixed_dims[i])
end


nc = mymdbm.ncubes[1]
corner = nc.corner
axes = mymdbm.axes
N = 3
Nf = 1
Nsubstep = 0


#
fixed_dims = Vector{Tuple{Int,Bool}}(undef, 0)
nc_template = mymdbm.T01

posall_tree = interpsubcubesolution([nc_template], [fixed_dims], corner, axes, Nf, Nc, N)
#faces, fixed_dims = generate_sub_faces(nc_template, fixed_dims)
#posall_tree = interpsubcubesolution(faces, fixed_dims, corner, axes, Nf, Nc, N)

function plot_tree(posall_tree,plotfun,argsDict,fighanle)
    plotfun(fighanle,( [i] for i in posall_tree[1][1] )...;plotargs...)
    for loc_plot in posall_tree[1][2] 
        plot_tree(loc_plot,plotfun,argsDict,fighanle)
    end
end


plotargs= Dict(:ms => 2)
fighanle=scatter()
plot_tree(posall_tree,scatter!,plotargs,fighanle)
scatter!()



