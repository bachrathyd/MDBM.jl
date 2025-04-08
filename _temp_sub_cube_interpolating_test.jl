using StaticArrays
using LinearAlgebra

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
- `fixed_dims::Vector{Int}`: Dimensions already fixed in the current face.

# Returns
- `sub_faces`: Vector containing sub-faces derived from the current face.
- `sub_fixed_dims`: Vector containing tuples indicating the new dimension and side fixed for each sub-face.
"""
function generate_sub_faces(face, fixed_dims::Vector{Int})
    n = length(face[1])
    sub_faces = []
    sub_fixed_dims = []
    free_dims = setdiff(1:n, fixed_dims)

    for dim in free_dims
        for side in [false, true]
            sub_face = filter(c -> c[dim] == side, face)
            push!(sub_faces, sub_face)
            push!(sub_fixed_dims, (dim, side))
        end
    end

    return sub_faces, sub_fixed_dims
end


# Example usage for a 3-cube:

N=3
faces, face_fixed_dims = generate_faces_indices(N)

println("Original faces:")
for (i, face) in enumerate(faces)
    println("Face $(i), fixed dim: ", face_fixed_dims[i])
    for corner in face
        println(" ", Int.(corner))
    end
end

N1 = 2
println("\nSub-faces of first face with fixed dimension [$N1]:")
sub_faces, sub_fixed_dims = generate_sub_faces(faces[N1], [face_fixed_dims[N1][1]])
for (i, sub_face) in enumerate(sub_faces)
    println("Sub-face $(i), fixed dim: ", sub_fixed_dims[i])
    for corner in sub_face
        println(" ", Int.(corner))
    end
end

N2 = 2
aa, bb = generate_sub_faces(sub_faces[N2], [face_fixed_dims[N1][1], sub_fixed_dims[N2][1]])

aa,bb=mymdbm.T01, zeros(Int,0)
#cc=[generate_sub_faces(aaloc,[bbloc[1]]) for (aaloc,bbloc) in zip(aa,bb)]

T11pinv = MDBM.T11pinvmaker(Val(N-1))
face = faces[1]

scatter(x_sol,y_sol,z_sol,ms=1)

#scatter()

nc=mymdbm.ncubes[1]
for nc in mymdbm.ncubes
    x0,y0,z0=getinterpolatedsolution(nc.posinterp, nc.corner, mymdbm.axes)
   # scatter!([x0],[y0],[z0],ms=5)
    #x0,y0=getinterpolatedsolution(nc.posinterp, nc.corner, mymdbm.axes)
    #scatter!([x0],[y0],ms=5)

   #[plot!(-1 .+ 2 .* getindex.(f, 1), -1 .+ 2 .* getindex.(f, 2), -1 .+ 2 .* getindex.(f, 3)) for f in faces]
   # scatter!()
   (face, fixdim) = (faces[1], face_fixed_dims[1])
    for (face, fixdim) in zip(faces, face_fixed_dims)
        #begin  (face, fixdim) = (faces[1], face_fixed_dims[1])
        posinterp = zeros(Float64, N)
        free_dims = setdiff(1:N, [fixdim[1]])
        posinterp[fixdim[1]] = -1 + 2 * fixdim[2]
        FunTupleVector = MDBM.getcornerval([nc.corner .+ nc.size .* T for T in face], mymdbm)
        #posinterp[free_dims] .= MDBM.fit_hyperplane(FunTupleVector, N - 1, 2, 0, typeof(FunTupleVector[1][1][1]), T11pinv)
        posinterp[free_dims] .= MDBM.fit_hyperplane(FunTupleVector, N - 1, 1, 0, typeof(FunTupleVector[1][1][1]), T11pinv)

        normp = 50000.0
        ncubetolerance = 0.001
        # @show posinterp
        # @show norm(posinterp, normp) < 1.0 + ncubetolerance
        if norm(posinterp, normp) < 1.0 + ncubetolerance
            x,y,z=getinterpolatedsolution(posinterp, nc.corner, mymdbm.axes)
            plot!([x0,x],[y0,y],[z0,z],lw=4)
            #x,y=getinterpolatedsolution(posinterp, nc.corner, mymdbm.axes)
            #scatter!([x],[y],[z])
            #plot!([x0,x],[y0,y],lw=4)
        end
    end
end
scatter!()


