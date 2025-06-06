function axdoubling!(ax::Axis)
    #WRONG: sort!(Base.append!(ax.ticks, ax.ticks[1:end - 1] + diff(ax.ticks) / 2); alg = QuickSort)
    #the axis could be in a non-sorted form

    doubleindex = (2:2*length(ax.ticks)) .÷ 2
    doubleindex[2:2:end] .+= length(ax.ticks)
    Base.append!(ax.ticks, ax.ticks[1:end-1] + diff(ax.ticks) / 2)
    ax.ticks[1:end] .= ax.ticks[doubleindex]
end


"""
    axesextend!(mdbm::MDBM_Problem{fcT,N,Nf,Nc,t01T,t11T,IT,FT,aT}, axisnumber::Integer; prepend::AbstractArray=[], append::AbstractArray=[])

Extend parameter space of the selected Axis by prepend and append the grid.

Note, that the resolution and the monotonically is not checked.

# Examples
extending the second parameter of the mdbm problem:
```jldoctest
julia> axesextend!(mymdbm,2,prepend=-6.2:0.05:-2.2,append=2.2:0.2:3);
```
"""
function axesextend!(mdbm::MDBM_Problem{fcT,N,Nf,Nc,t01T,t11T,IT,FT,aT}, axisnumber::Integer; prepend::AbstractArray=[], append::AbstractArray=[]) where {fcT,N,Nf,Nc,t01T,t11T,IT,FT,aT}
    preplength = length(prepend)
    if preplength > 0
        prepend!(mdbm.axes[axisnumber].ticks, prepend)
        @inbounds for nc in mdbm.ncubes
            nc.corner[axisnumber] += preplength
        end
    end
    if length(append) > 0
        Base.append!(mdbm.axes[axisnumber].ticks, append)
    end
end

# function corner(nc::NCube{IT,FT,N,Nfc}, T01)::Vector{MVector} where IT where FT where N
function corner(nc::NCube{IT,FT,N,Nfc}, T01) where {IT,FT,N,Nfc}
    [nc.corner .+ nc.size .* T for T in T01]
end


# function corner(ncubes::Vector{NCube{IT,FT,N,Nfc}}, T01)::Vector{Vector{MVector}} where IT where FT where N
function corner(ncubes::Vector{<:NCube}, T01) #{IT,FT,N,Nfc} where {IT,FT,N,Nfc}
    [corner(nc, T01) for nc in ncubes]
    # [nc.corner .+ nc.size .* T for nc in ncubes for T in T01]
end

# function corner(mdbm::MDBM_Problem{fcT,N,Nf,Nc,t01T,t11T,IT,FT,aT})::Vector{Vector{SVector}}  where fcT where N where Nf where Nc
function corner(mdbm::MDBM_Problem{fcT,N,Nf,Nc,t01T,t11T,IT,FT,aT}) where {fcT,N,Nf,Nc,t01T,t11T,IT,FT,aT}
    [corner(nc, mdbm.T01) for nc in mdbm.ncubes]
end



function getcornerval(mdbm::MDBM_Problem{fcT,N,Nf,Nc,t01T,t11T,IT,FT,aT}) where {fcT,N,Nf,Nc,t01T,t11T,IT,FT,aT}#get it for all
    #getcornerval.(mdbm.ncubes,Ref(mdbm))

    funargs = map(x -> ((mdbm.axes))(x...), Base.Iterators.flatten(corner(mdbm)))
    #mdbm.fc(unique(funargs))#prcomputed if it was not done before
    mdbm.fc.(funargs)

    #map(x -> ((mdbm.fc) ∘ (mdbm.axes))(x...), Base.Iterators.flatten(corner(mdbm)))
end

function getcornerval(ncube::NCube{IT,FT,N,Nfc}, mdbm::MDBM_Problem{fcT,N,Nf,Nc,t01T,t11T,IT,FT,aT}) where {fcT,N,Nf,Nc,t01T,t11T,IT,FT,aT,Nfc}
    #PostionsVectors=(mdbm.axes.(corner(nc,mdbm.T01)))
    # (mdbm.fc).( (mdbm.axes).(corner(nc,mdbm.T01)) )

    funargs = map(x -> ((mdbm.axes))(x...), corner(ncube, mdbm.T01))
    #mdbm.fc(unique(funargs))#prcomputed if it was not done before
    mdbm.fc.(funargs)

    # map(x -> ((mdbm.fc) ∘ (mdbm.axes))(x...), corner(ncubes, mdbm.T01))
end

function getcornerval(corners, mdbm::MDBM_Problem{fcT,N,Nf,Nc,t01T,t11T,IT,FT,aT}) where {fcT,N,Nf,Nc,t01T,t11T,IT,FT,aT}
    funargs = map(x -> ((mdbm.axes))(x...), corners)
    #mdbm.fc(unique(funargs))#prcomputed if it was not done before
    mdbm.fc.(funargs)

    # map(x -> ((mdbm.fc) ∘ (mdbm.axes))(x...), corners)
end


function T01maker(valk::Val{kdim}) where {kdim}
    SVector{2^kdim}([
        SVector{kdim}([isodd(x ÷ (2^y)) for y in 0:(kdim-1)])
        for x in 0:(2^kdim-1)])
end

function T11pinvmaker(valk::Val{kdim}) where {kdim}
    T11pinv = ([isodd(x ÷ (2^y)) for y in 0:kdim, x in 0:(2^kdim-1)] * 2.0 .- 1.0) / (2^kdim)
    SMatrix{kdim + 1,2^kdim}(T11pinv)
end


# #TODO : deprecated - same results can be achieved with the generate_sub_faces function
# """
#     generate_faces_indices(n)
# 
# Generates the faces of an n-dimensional cube.
# 
# # Arguments
# - `n::Int`: Dimension of the cube.
# 
# # Returns
# - `faces`: Vector containing the faces, each face represented by a vector of corners.
# - `fixed_dims`: Vector containing tuples indicating the dimension and the side (0 or 1) fixed for each face.
# """
# function generate_faces_indices(n)
#     corners = [SVector{n}(digits(i, base=2, pad=n) .== 1) for i in 0:(2^n-1)]
#     faces = Vector{Vector{SVector{n,Bool}}}(undef, 2 * n)
#     fixed_dims = Vector{Tuple{Int,Bool}}(undef, 2 * n)
# 
#     idx = 1
#     for dim in 1:n
#         for side in [false, true]
#             face = filter(c -> c[dim] == side, corners)
#             faces[idx] = face
#             fixed_dims[idx] = (dim, side)
#             idx += 1
#         end
#     end
# 
#     return faces, fixed_dims
# end

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

    @inbounds for dim in free_dims
        for side in [false, true]
            sub_face = filter(c -> c[dim] == side, face)
            push!(sub_faces, sub_face)
            push!(sub_fixed_dims, (dim, side))
        end
    end

    all_fixed_dims = [vcat(fixed_dims, [sf]) for sf in sub_fixed_dims]

    return sub_faces, all_fixed_dims#, sub_fixed_dims, 
end



function interpsubcubesolution!(mdbm::MDBM_Problem{fcT,N,Nf,Nc,t01T,t11T,IT,FT,aT}; normp=20.0, ncubetolerance=0.5) where {fcT,N,Nf,Nc,t01T,t11T,IT,FT,aT}
    for nc in mdbm.ncubes
        interpsubcubesolution!(nc, mdbm, normp=normp, ncubetolerance=ncubetolerance)
    end
    return nothing
end
function interpsubcubesolution!(ncubes::Vector{<:NCube}, mdbm::MDBM_Problem{fcT,N,Nf,Nc,t01T,t11T,IT,FT,aT}; normp=20.0, ncubetolerance=0.5) where {fcT,N,Nf,Nc,t01T,t11T,IT,FT,aT}
    for nc in ncubes
        interpsubcubesolution!(nc, mdbm, normp=normp, ncubetolerance=ncubetolerance)
    end
    return nothing
end

function interpsubcubesolution!(nc::NCube, mdbm::MDBM_Problem{fcT,N,Nf,Nc,t01T,t11T,IT,FT,aT}; normp=20.0, ncubetolerance=0.5) where {fcT,N,Nf,Nc,t01T,t11T,IT,FT,aT}

    fixed_dims = Vector{Tuple{Int,Bool}}(undef, 0)
    nc_template = mdbm.T01

    # posall_tree = MDBM.PositionTree(nc.posinterp)
    edges, all_edge_fixed_dims = MDBM.generate_sub_faces(nc_template, fixed_dims)

    MDBM.interpsubcubesolution!(nc.posinterp, edges, all_edge_fixed_dims, nc.corner, nc.size, mdbm, normp=normp, ncubetolerance=ncubetolerance)

end

function interpsubcubesolution!(posall_tree, faces, fixed_dims, corner, size, mdbm::MDBM_Problem{fcT,N,Nf,Nc,t01T,t11T,IT,FT,aT}; normp=20.0, ncubetolerance=0.5) where {fcT,N,Nf,Nc,t01T,t11T,IT,FT,aT}


    @inbounds for (face, fixdim) in zip(faces, fixed_dims)
        FunTupleVector = MDBM.getcornerval([corner .+ size .* T for T in face], mdbm)

        fixed_dimensions_only = [fd[1] for fd in fixdim]
        free_dims = setdiff(1:N, fixed_dimensions_only)
        Nfree = length(free_dims)

        T11pinv = MDBM.T11pinvmaker(Val(length(free_dims)))

        posinterp = zeros(MVector{N,FT})
        grad = zeros(MMatrix{N,Nf + Nc,FT})

        # As = MVector{Nf + Nc,FT}(undef)
        # ns = MMatrix{length(free_dims),Nf + Nc,FT}(undef)
        # fit_hyperplane!(FunTupleVector, Val(Nfree), Val(Nf), Val(Nc), T11pinv,
        # posinterp[free_dims],grad[free_dims,:],As,ns)

        p, g = fit_hyperplane(FunTupleVector, Val(Nfree), Val(Nf), Val(Nc), typeof(FunTupleVector[1][1][1]), T11pinv)
        posinterp[free_dims] .= p

        posinterp[getindex.(fixdim, 1)] .= -1 .+ 2 .* getindex.(fixdim, 2)


        #normp = 50000.0
        #ncubetolerance = 0.2

        if norm(posinterp, normp) < 1.0 + ncubetolerance
            #print("$Nfree ok:")
            #@show face

            #push!(posall_tree.subpoints, PositionTree(getinterpolatedsolution(posinterp, corner, mdbm.axes)))
            push!(posall_tree.subpoints, PositionTree(posinterp))
            if Nfree > Nf
                edges, all_edge_fixed_dims = generate_sub_faces(face, fixdim)
                interpsubcubesolution!(posall_tree.subpoints[end], edges, all_edge_fixed_dims, corner, size, mdbm, normp=normp, ncubetolerance=ncubetolerance)
            end
        else
            #println("$Nfree x")
        end
    end
    return nothing
    # return posall_tree
end


function extract_paths_local(nc::NCube)
    extract_paths(nc.posinterp)
end
function extract_paths(nc::NCube, mdbm::MDBM_Problem{fcT,N,Nf,Nc,t01T,t11T,IT,FT,aT}) where {fcT,N,Nf,Nc,t01T,t11T,IT,FT,aT}
    [getinterpolatedsolution.(posloc, Ref(nc.corner), Ref(mdbm.axes)) for posloc in extract_paths_local(nc)]
end

function extract_paths(mdbm::MDBM_Problem{fcT,N,Nf,Nc,t01T,t11T,IT,FT,aT}) where {fcT,N,Nf,Nc,t01T,t11T,IT,FT,aT}
    extract_paths.(mdbm.ncubes, Ref(mdbm))
end



"""
    extract_paths(tree::PositionTree{N,T})

Extracts all root-to-leaf paths from a PositionTree.

# Arguments
- `tree::PositionTree{N,T}`: The root of the PositionTree.

# Returns
- An array of paths, each path being a vector of positions (as SVectors) from root to leaf.

# Example
```julia
root = PositionTree([1.0, 2.0, 3.0])
push!(root.subpoints, PositionTree([4.0, 5.0, 6.0]))
push!(root.subpoints, PositionTree([7.0, 8.0, 9.0]))

paths = extract_paths(root)
```
"""
function extract_paths(tree::PositionTree{N,T}, current_path=Vector{SVector{N,T}}()) where {N,T}
    paths_loc = Vector{Vector{SVector{N,T}}}()

    # Append current node's position to the path
    #new_path = [current_path; SVector{N,T}(tree.p)]
    new_path = deepcopy(current_path)
    push!(new_path, tree.p)

    # If the node is a leaf, return the accumulated path
    if isempty(tree.subpoints)
        push!(paths_loc, new_path)
        #println("Leaf node added")
    else
        # Otherwise, recursively explore subpoints

        #println("--------~~~~~~~~~~~~~~~~--------------")
        #println("Subtree length: ", length(tree.subpoints))
        for subtree in tree.subpoints
            #@show (subtree, new_path)
            Base.append!(paths_loc, extract_paths(subtree, new_path))
        end
    end

    return paths_loc
end

function issingchange(FunTupleVector::AbstractVector, Nf::Integer, Nc::Integer)::Bool
    @inbounds all([
        [
            any((c) -> !isless(c[1][fi], zero(c[1][fi])), FunTupleVector)
            for fi in 1:Nf#length(FunTupleVector[1][1])
        ]#any positive (or zeros) value (!isless)
        [
            any((c) -> !isless(zero(c[1][fi]), c[1][fi]), FunTupleVector)
            for fi in 1:Nf#length(FunTupleVector[1][1])
        ]#anany negative (or zero) value
        [
            any((c) -> !isless(c[2][fi], zero(c[2][fi])), FunTupleVector)
            for fi in 1:Nc#length(FunTupleVector[1][2])
        ]#any positive (or zeros) value (!isless)
    ])#check for all condition
end

function _interpolate!(ncubes::Vector{<:NCube}, mdbm::MDBM_Problem{fcT,N,Nf,Nc,t01T,t11T,IT,FT,aT}, ::Type{Val{0}}, normp=20.0, ncubetolerance=0.5, doThreadprecomp=true) where {fcT,N,Nf,Nc,t01T,t11T,IT,FT,aT}
    # println("precompute the missing elemnts toghether!")
    funargs = map(x -> ((mdbm.axes))(x...), Base.Iterators.flatten(corner(mdbm)))
    if doThreadprecomp && typeof(mdbm.fc) <: MemF  #TODO: It is parallel, but Somehow it is slower!!
        mdbm.fc(unique!(funargs))#prcomputed if it was not done before
    else
        # mdbm.fc.(unique(funargs))#prcomputed if it was not done before
    end
    isbracketing = map(nc -> issingchange(getcornerval(nc, mdbm), Nf, Nc), ncubes) #TODO: parallelize
    deleteat!(ncubes, .!isbracketing)#removeing the non-bracketing ncubes

    @inbounds for nc in ncubes
        nc.posinterp.p[:] .= zero(FT)
    end
    return nothing
end

function fit_hyperplane_CHATGPT!(
    FunTupleVector,
    ::Val{Nfree}, ::Val{Nf}, ::Val{Nc}, T11pinv,
    posinterp::AbstractVector{FT},
    grad::AbstractMatrix{FT},
) where {Nfree,Nf,Nc,FT}

    # number of corners = 2^Nfree
    #const 
    C = 1 << Nfree

    # reuse these stack‐allocated buffers:
    vals = MVector{C,FT}(undef)                  # to hold corner‐values
    A_sf = MMatrix{C,C,FT}(undef)              # temporary, if you need it
    # but we’re going to go straight to a ns x (Nf+Nc) matrix:
    ns_sf = MMatrix{Nfree,Nf + Nc,FT}(undef)     # gradients
    As_sf = MVector{Nf + Nc,FT}(undef)             # intercepts

    # 1) fill As_sf and ns_sf via static loops
    # -----------------------------------------------------
    # for each equation fᵢ or active constraint cⱼ,
    # build the C-vector of values, multiply by T11pinv once,
    # split off last entry into As_sf, the rest into ns_sf:
    eq = 0
    @inbounds for k in 1:Nf
        @inbounds for x in 1:C
            vals[x] = FunTupleVector[x][1][k]
        end
        let sol = T11pinv * vals  # this is all static
            As_sf[k] = sol[C]
            for d in 1:Nfree
                ns_sf[d, k] = sol[d]
            end
        end
        eq += 1
    end

    if Nc > 0
        @inbounds for k in 1:Nc
            # only if this constraint is “active” in this cube
            if any(c -> !isless(zero(c[2][k]), c[2][k]), FunTupleVector)
                @inbounds for x in 1:C
                    vals[x] = FunTupleVector[x][2][k]
                end
                let sol = T11pinv * vals
                    As_sf[Nf+eq] = sol[C]
                    for d in 1:Nfree
                        ns_sf[d, Nf+eq] = sol[d]
                    end
                end
                eq += 1
            end
        end
    end

    # 2) solve the small normal‐equation system via a single static \
    # -----------------------------------------------------
    # a is Nfree×eq,  d is eq-vector:
    a = @view ns_sf[:, 1:eq]
    d = view(As_sf, 1:eq)

    # use back-slash on a static array (gives you the least-squares solution)
    # it will pick QR or Cholesky under the hood, with _no_ SVD allocation
    try
        p = a \ d   # MVector{Nfree,FT}
    catch
        # singular → push it out of the cube
        fill!(posinterp, 1e3)
        fill!(grad, NaN)
        return
    end

    # 3) copy back into the user buffers and return
    @inbounds for i in 1:Nfree
        posinterp[i] = p[i]
        for j in 1:eq
            grad[i, j] = ns_sf[i, j]
        end
    end

    return
end


function fit_hyperplane!(FunTupleVector, ::Val{N}, ::Val{Nf}, ::Val{Nc}, T11pinv,
    posinterp::AbstractVector{FT},
    grad::AbstractMatrix{FT},
    As::AbstractVector{FT}, ns::AbstractMatrix{FT}) where {N,Nf,Nc,FT}
    # number of corners = 2^Nfree


    #TODO: store them somewhere - to reduce memory allcoations
    #As = MVector{Nf + Nc,FT}(undef)
    #ns = MMatrix{N,Nf + Nc,FT}(undef)
    @inbounds if Nc == 0 || all(
        any((c) -> !isless(c[2][fi], zero(c[2][fi])), FunTupleVector)
        for fi in 1:length(FunTupleVector[1][2])
    )# check the constraint: do wh have to compute at all?!?


        #for f---------------------
        for kf = 1:Nf#length(FunTupleVector[1][1])
            solloc = T11pinv * [FunTupleVector[kcubecorner][1][kf] for kcubecorner = 1:length(FunTupleVector)]
            ns[:, kf] .= solloc[1:end-1]
            As[kf] = solloc[end]

        end

        #for c---if needed: that is, it is close to the boundary--------
        activeCostraint = 0
        for kf = 1:Nc#length(FunTupleVector[1][2])
            #TODO: mivan a többszörös C teljesülése esetén!?!??!# ISSUE
            if any((c) -> !isless(zero(c[2][kf]), c[2][kf]), FunTupleVector) && length(As) <= N  # use a constraint till it reduces the dimension to zero (point) and no further
                solloc = T11pinv * [FunTupleVector[kcubecorner][2][kf] for kcubecorner = 1:length(FunTupleVector)]
                activeCostraint += 1
                As[Nf+activeCostraint] = solloc[end]
                ns[:, Nf+activeCostraint] .= solloc[1:end-1]
            end
        end

        # TOO sepcial, newer should be happen
        #if (Nf + activeCostraint) == 0 
        #    posinterp = zeros(MVector{N,FT})
        # first two leads to error in the fucntion of constant due to the behaviour of pinv (eg. [1 0;0 1;0 0]/[1 0; 0 0; 0 0])
        # nc.posinterp[:] .= transpose(view(ns, :, 1:(Nf + activeCostraint))) \ view(As, 1:(Nf + activeCostraint));#TEST
        # nc.posinterp[:] .= transpose(ns[:, 1:(Nf + activeCostraint)]) \ As[1:(Nf + activeCostraint)];#TEST
        # nc.posinterp[:] .= ns[:, 1:(Nf + activeCostraint)] * ((transpose(ns[:, 1:(Nf + activeCostraint)]) * ns[:, 1:(Nf + activeCostraint)]) \ view(As, 1:(Nf + activeCostraint)));

        a = @view ns[:, 1:(Nf+activeCostraint)]
        d = view(As, 1:(Nf+activeCostraint))
        A = transpose(a) * a
        if rank(A) == minimum(size(A))
            posinterp[:] .= a * (A \ d)
        else
            posinterp[:] .= 1000.0#put it outside the cube!
        end

        #end
        grad = ns
        #for c---------------------
    else
        posinterp[:] .= 1000.0#put it outside the cube!
        grad[:] .= NaN#put it outside the cube!
    end
    return posinterp, grad
end

function fit_hyperplane(FunTupleVector, ::Val{N}, ::Val{Nf}, ::Val{Nc}, FT, T11pinv) where {N,Nf,Nc}
    # number of corners = 2^Nfree

    posinterp = zeros(MVector{N,FT})
    grad = zeros(MMatrix{N,Nf + Nc,FT})

    #TODO: store them somewhere - to reduce memory allcoations
    As = MVector{Nf + Nc,FT}(undef)
    ns = MMatrix{N,Nf + Nc,FT}(undef)
    @inbounds if Nc == 0 || all(
        any((c) -> !isless(c[2][fi], zero(c[2][fi])), FunTupleVector)
        for fi in 1:length(FunTupleVector[1][2])
    )# check the constraint: do wh have to compute at all?!?


        #for f---------------------
        for kf = 1:Nf#length(FunTupleVector[1][1])
            solloc = T11pinv * [FunTupleVector[kcubecorner][1][kf] for kcubecorner = 1:length(FunTupleVector)]
            As[kf] = solloc[end]
            ns[:, kf] .= solloc[1:end-1]
        end

        #for c---if needed: that is, it is close to the boundary--------
        activeCostraint = 0
        for kf = 1:Nc#length(FunTupleVector[1][2])
            #TODO: mivan a többszörös C teljesülése esetén!?!??!# ISSUE
            if any((c) -> !isless(zero(c[2][kf]), c[2][kf]), FunTupleVector) && length(As) <= N  # use a constraint till it reduces the dimension to zero (point) and no further
                solloc = T11pinv * [FunTupleVector[kcubecorner][2][kf] for kcubecorner = 1:length(FunTupleVector)]
                activeCostraint += 1
                As[Nf+activeCostraint] = solloc[end]
                ns[:, Nf+activeCostraint] .= solloc[1:end-1]
            end
        end

        if (Nf + activeCostraint) == 0
            posinterp = zeros(MVector{N,FT})
        else
            # first two leads to error in the fucntion of constant due to the behaviour of pinv (eg. [1 0;0 1;0 0]/[1 0; 0 0; 0 0])
            # nc.posinterp[:] .= transpose(view(ns, :, 1:(Nf + activeCostraint))) \ view(As, 1:(Nf + activeCostraint));#TEST
            # nc.posinterp[:] .= transpose(ns[:, 1:(Nf + activeCostraint)]) \ As[1:(Nf + activeCostraint)];#TEST
            # nc.posinterp[:] .= ns[:, 1:(Nf + activeCostraint)] * ((transpose(ns[:, 1:(Nf + activeCostraint)]) * ns[:, 1:(Nf + activeCostraint)]) \ view(As, 1:(Nf + activeCostraint)));

            a = ns[:, 1:(Nf+activeCostraint)]
            d = view(As, 1:(Nf+activeCostraint))
            A = transpose(a) * a
            if rank(A) == minimum(size(A))
                posinterp[:] .= a * (A \ d)
            else
                posinterp[:] .= 1000.0#put it outside the cube!
            end

        end
        grad = ns
        #for c---------------------
    else
        posinterp[:] .= 1000.0#put it outside the cube!
        grad[:] .= NaN#put it outside the cube!
    end
    return posinterp, grad
end

function _interpolate!(ncubes::Vector{<:NCube}, mdbm::MDBM_Problem{fcT,N,Nf,Nc,t01T,t11T,IT,FT,aT}, ::Type{Val{1}}; normp=20.0, ncubetolerance=0.5, doThreadprecomp=true) where {fcT,N,Nf,Nc,t01T,t11T,IT,FT,aT}
    # println("precompute the missing elemnts toghether!")
    funargs = map(x -> ((mdbm.axes))(x...), Base.Iterators.flatten(corner(ncubes, mdbm.T01)))
    if doThreadprecomp && typeof(mdbm.fc) <: MemF #Memoization is used #TODO: It is parallel, but Somehow it is slower!!
        #mdbm.fc(unique(funargs))#prcomputed if it was not done before
        mdbm.fc(unique!(funargs))#prcomputed if it was not done before
    else
        mdbm.fc.(unique!(funargs))#prcomputed if it was not done before -> necessary for the later threaded for loop
    end

    #As = MVector{Nf + Nc,FT}(undef)
    #ns = MMatrix{N,Nf + Nc,FT}(undef)
    #for nc in ncubes
    Threads.@threads for i_nc in 1:length(ncubes)
        nc = ncubes[i_nc]
        As = MVector{Nf + Nc,FT}(undef)
        ns = MMatrix{N,Nf + Nc,FT}(undef)
        FunTupleVector = getcornerval(nc, mdbm)
        fit_hyperplane!(FunTupleVector, Val(N), Val(Nf), Val(Nc), mdbm.T11pinv,
            nc.posinterp.p, nc.gradient, As, ns)
        # nc.posinterp.p[:], nc.gradient[:] = fit_hyperplane(FunTupleVector, Val(N), Val(Nf), Val(Nc), FT, mdbm.T11pinv)

    end

    #TODO: what if it falls outside of the n-cube, it should be removed ->what shall I do with the bracketing cubes?
    # Let the user define it
    # filter!((nc)->sum((abs.(nc.posinterp)).^10.0)<(1.5 ^10.0),mdbm.ncubes)#1e6~=(2.0 ^20.0)


    filter!((nc) -> norm(nc.posinterp.p, normp) < 1.0 + ncubetolerance, ncubes)
    #filter!((nc)->!any(isnan.(nc.posinterp.p)),mdbm.ncubes)

    return nothing
end

function _interpolate!(ncubes::Vector{<:NCube}, mdbm::MDBM_Problem{fcT,N,Nf,Nc,t01T,t11T,IT,FT,aT}, ::Type{Val{Ninterp}}; normp=20.0, ncubetolerance=0.5, doThreadprecomp=true) where {fcT,IT,FT,N,Ninterp,Nf,Nc,t01T,t11T,aT}
    error("order $(Ninterp) interpolation is not supperted (yet)")
end

"""
    interpolate!(mdbm::MDBM_Problem{fcT,N,Nf,Nc,t01T,t11T,IT,FT,aT}; interpolationorder::Int=1)

Interpolate the solution within the n-cube with order `interpolationorder`.
 - `interpolationorder=0` provide the midpoint of the n-cube
 - `interpolationorder=1` preform a linear fit and provide closest solution point to the mindpoint of the n-cube
"""
function interpolate!(mdbm::MDBM_Problem{fcT,N,Nf,Nc,t01T,t11T,IT,FT,aT}; interpolationorder::Int=1, normp=20.0, ncubetolerance=0.5, doThreadprecomp=true) where {fcT,N,Nf,Nc,t01T,t11T,IT,FT,aT}
    _interpolate!(mdbm.ncubes, mdbm, Val{interpolationorder}, normp=normp, ncubetolerance=ncubetolerance, doThreadprecomp=doThreadprecomp)
end



"""
    refine!(mdbm::MDBM_Problem{fcT,N,Nf,Nc,t01T,t11T,IT,FT,aT}; directions::Vector{T}=collect(Int64,1:N)) where N where Nf where Nc where T<:Integer

Double the resolution of the axes along the provided `directions` then halving the `ncubes` accordingly.

# Examples
```jldoctest
julia> refine!(mymdbm)
julia> refine!(mymdbm,[1,1,2])
```
"""
function refine!(mdbm::MDBM_Problem{fcT,N,Nf,Nc,t01T,t11T,IT,FT,aT}; directions::Vector{T}=collect(Int64, 1:N)) where {T<:Integer,fcT,IT,FT,N,Nf,Nc,t01T,t11T,aT}
    doubling!(mdbm, directions)
    refinencubes!(mdbm.ncubes, directions)
    return nothing
end

"""
    doubling!(mdbm::MDBM_Problem{fcT,N,Nf,Nc,t01T,t11T,IT,FT,aT}, directions::Vector{T}) where T<:Integer  where N where Nf where Nc

Double the resolution of the axes along the provided `directions`, then set the new size of the n-cubes.

# Examples
```jldoctest
julia> doubling!(mymdbm)
julia> doubling!(mymdbm,[1,2])
```
"""
function doubling!(mdbm::MDBM_Problem{fcT,N,Nf,Nc,t01T,t11T,IT,FT,aT}, directions::Vector{T}) where {T<:Integer} where {fcT,N,Nf,Nc,t01T,t11T,IT,FT,aT}
    axdoubling!.(mdbm.axes[directions])
    @inbounds for nc in mdbm.ncubes
        for dir in directions
            nc.corner[dir] = (nc.corner[dir] - 1) * 2 + 1
            nc.size[dir] *= 2
        end
    end
end

"""
     refinencubes!(ncubes::Vector{<:NCube}, directions::Vector{T})

Halving the `ncubes` along the provided `directions`.

# Examples
```jldoctest
julia> refinencubes!(mymdbm.ncubes)
julia> refinencubes!(mymdbm.ncubes,mymdbm,[1,2])
```
"""
# {IT,FT,N,Nfc} where {IT,FT,N,Nfc} 
function refinencubes!(ncubes::Vector{<:NCube{IT,FT,N,Nfc}}, directions::Vector{T}) where {T<:Integer,IT,FT,N,Nfc} #where t01T # where IT where FT where N
    sizehint!(ncubes, (length(ncubes) * 2^length(directions)) * 1.2) # preallocate memory for the new ncubes, and some overhead, to have extra space fo neighbour chacking
    @inbounds for dir in directions
        for nc in ncubes
            nc.size[dir] /= 2
        end
        NumofNCubes = length(ncubes)
        #        Base.append!(ncubes, deepcopy(ncubes))
        #        for nc in ncubes[1:NumofNCubes]
        #            nc.corner[dir] = nc.corner[dir] + nc.size[dir]
        #        end

        # Pre‐allocate a Vector to hold the new cubes:
        resize!(ncubes, 2 * NumofNCubes)         # now capacity is 2*old_n


        #new_cubes = Vector{NCube{IT,FT,N,Nfc}}(undef, NumofNCubes)
        old = deepcopy(ncubes[1])
        for i in 1:NumofNCubes
            old = ncubes[i]
            # 1) make a fresh corner vector, shifted by old.size[dir]:
            corner2 = MVector{N,IT}(old.corner)      # copy of old.corner
            corner2[dir] += old.size[dir]
            # 2) make a fresh size vector (you can just copy it):
            size2 = MVector{N,IT}(old.size)
            # 3) initialize a brand‐new PositionTree (all zeros) and gradient matrix:
            pos2 = PositionTree(zeros(MVector{N,FT}))        # or whatever default you need
            grad2 = MMatrix{N,Nfc,FT}(undef)
            # 4) construct the new cube with the same bracketing flag (or set it if desired):
            ncubes[i+NumofNCubes] = NCube{IT,FT,N,Nfc}(corner2, size2, pos2, true, grad2)
        end
        # finally, append all the new ones at once:

        #new_all_cubes=merge_sorted(ncubes, new_cubes)
        #empty!(ncubes)
        #append!(ncubes, new_all_cubes) # no need for sort!-but slower, do to the memory allcoation

        #append!(ncubes, new_cubes)

        half_merge_with_buffer!(ncubes)

    end
    #sort!(ncubes; alg=QuickSort)
    return nothing
end



"""
    getinterpolatedsolution(mdbm::MDBM_Problem{fcT,N,Nf,Nc,t01T,t11T,IT,FT,aT})

Provide the interpolated coordinates of the all the detected solution (approximately where foo(x,y) == 0 and c(x,y)>0).
"""
function getinterpolatedsolution(mdbm::MDBM_Problem{fcT,N,Nf,Nc,t01T,t11T,IT,FT,aT}) where {fcT,N,Nf,Nc,t01T,t11T,IT,FT,aT}
    getinterpolatedsolution(mdbm.ncubes, mdbm)
end

"""
    getinterpolatedsolution(ncubes::Vector{<:NCube},mdbm::MDBM_Problem{N,Nf,Nc})

Provide the interpolated coordinates of the detected solution for the selected n-cubes (approximately where foo(x,y) == 0 and c(x,y)>0).
"""
function getinterpolatedsolution(ncubes::Vector{<:NCube}, mdbm::MDBM_Problem{fcT,N,Nf,Nc,t01T,t11T,IT,FT,aT}) where {fcT,N,Nf,Nc,t01T,t11T,IT,FT,aT}
    @inbounds [
        [
            (typeof(mdbm.axes[i].ticks).parameters[1])(
                (mdbm.axes[i].ticks[nc.corner[i]] * (1.0 - (nc.posinterp.p[i] + 1.0) / 2.0) +
                 mdbm.axes[i].ticks[nc.corner[i]+nc.size[i]] * ((nc.posinterp.p[i] + 1.0) / 2.0))
            )
            for nc in ncubes]#::Vector{typeof(mdbm.axes[i].ticks).parameters[1]}
        for i in 1:length(mdbm.axes)]
end

"""
    getinterpolatedsolution(nc::NCube{IT,FT,N,Nfc}, mdbm::MDBM_Problem{fcT,N,Nf,Nc,t01T,t11T,IT,FT,aT})

Provide the interpolated coordinates of the solution inside the provided n-cube `nc` (approximately where foo(x,y) == 0 and c(x,y)>0).
"""
function getinterpolatedsolution(nc::NCube{IT,FT,N,Nfc}, mdbm::MDBM_Problem{fcT,N,Nf,Nc,t01T,t11T,IT,FT,aT}) where {fcT,N,Nf,Nc,t01T,t11T,IT,FT,aT,Nfc}
    @inbounds [
        (typeof(mdbm.axes[i].ticks).parameters[1])(
            (mdbm.axes[i].ticks[nc.corner[i]] * (1.0 - (nc.posinterp.p[i] + 1.0) / 2.0) +
             mdbm.axes[i].ticks[nc.corner[i]+nc.size[i]] * ((nc.posinterp.p[i] + 1.0) / 2.0)))#::Vector{typeof(mdbm.axes[i].ticks).parameters[1]}
        for i in 1:length(mdbm.axes)]
end

function getinterpolatedsolution(posinterp, corner, axes::Axes{N,aT}) where {N,aT}
    @inbounds [
        (typeof(axes[i].ticks).parameters[1])((axes[i].ticks[corner[i]] * (1.0 - (posinterp[i] + 1.0) / 2.0) +
                                               axes[i].ticks[corner[i]+1] * ((posinterp[i] + 1.0) / 2.0)))#::Vector{typeof(mdbm.axes[i].ticks).parameters[1]}
        for i in 1:length(axes)]
end


function getinterpolatedgradient(ncubes::Vector{<:NCube}, mdbm::MDBM_Problem{fcT,N,Nf,Nc,t01T,t11T,IT,FT,aT}) where {fcT,N,Nf,Nc,t01T,t11T,IT,FT,aT}
    #println("getinterpolatedgradient")
    @inbounds [[
        [
            nc.gradient[i, fi] / ((mdbm.axes[i].ticks[nc.corner[i]+nc.size[i]] -
                                   mdbm.axes[i].ticks[nc.corner[i]]) / 2.0)
            for nc in ncubes]#::Vector{typeof(mdbm.axes[i].ticks).parameters[1]}
        for i in 1:length(mdbm.axes)] for fi in 1:size(ncubes[1].gradient, 2)]
end


function getinterpolatedgradient(nc::NCube, mdbm::MDBM_Problem{fcT,N,Nf,Nc,t01T,t11T,IT,FT,aT}) where {fcT,N,Nf,Nc,t01T,t11T,IT,FT,aT}
    @inbounds [[
        (mdbm.axes[i].ticks[nc.corner[i]] * (1.0 .- (nc.gradient[i, fi] .+ 1.0) ./ 2.0) +
         mdbm.axes[i].ticks[nc.corner[i]+nc.size[i]] * ((nc.gradient[i, fi] .+ 1.0) ./ 2.0))
        for i in 1:length(mdbm.axes)] for fi in 1:size(nc.gradient, 2)]
end

"""
    getevaluatedpoints(mdbm::MDBM_Problem{fcT,N,Nf,Nc,t01T,t11T,IT,FT,aT})

Provide all the coordinates where the function is evaluated.
"""
function getevaluatedpoints(mdbm::MDBM_Problem{fcT,N,Nf,Nc,t01T,t11T,IT,FT,aT}) where {fcT,N,Nf,Nc,t01T,t11T,IT,FT,aT}
    #    [[x.callargs[i] for x in mdbm.fc.fvalarg] for i in 1:N]
    [getindex.(mdbm.fc.fvalarg.keys, i) for i in 1:N]
end

"""
    getevaluatedfunctionvalues(mdbm::MDBM_Problem{fcT,N,Nf,Nc,t01T,t11T,IT,FT,aT})

Provide the functionvalues for all the coordinates where the function is evaluated.
(see: `getevaluatedfunctionvalues`)
"""
function getevaluatedfunctionvalues(mdbm::MDBM_Problem{fcT,N,Nf,Nc,t01T,t11T,IT,FT,aT}) where {fcT,N,Nf,Nc,t01T,t11T,IT,FT,aT}
    #[x.funval for x in mdbm.fc.fvalarg]
    getindex.(getindex.(Ref(mdbm.fc.fvalarg), mdbm.fc.fvalarg.keys), 1)
end

"""
    getevaluatedconstraintvalues(mdbm::MDBM_Problem{fcT,N,Nf,Nc,t01T,t11T,IT,FT,aT})

Provide the constraints values for all the coordinates where the function is evaluated.
(see: `getevaluatedconstraintvalues`)
"""
function getevaluatedconstraintvalues(mdbm::MDBM_Problem{fcT,N,Nf,Nc,t01T,t11T,IT,FT,aT}) where {fcT,N,Nf,Nc,t01T,t11T,IT,FT,aT}
    #[x.cval for x in mdbm.fc.fvalarg]
    getindex.(getindex.(Ref(mdbm.fc.fvalarg), mdbm.fc.fvalarg.keys), 2)
end

# function generateneighbours(ncubes::Vector{<:NCube},mdbm::MDBM_Problem{fcT,N,Nf,Nc,t01T,t11T,IT,FT,aT}) where N  where Nf where Nc
function generateneighbours(ncubes::Vector{<:NCube}, mdbm::MDBM_Problem{fcT,N,Nf,Nc,t01T,t11T,IT,FT,aT}) where {fcT,N,Nf,Nc,t01T,t11T,IT,FT,aT}
    #  #TODO: Let the user select the method
    # #-------all faceing/cornering neightbours----------------
    #  #neighbourind=1:2^length(mdbm.axes) #cornering neightbours also - unnecessary
    #  neighbourind=2 .^(0:(length(mdbm.axes)-1)) .+ 1 #neighbour only on the side
    #  T101=[-mdbm.T01[neighbourind]...,mdbm.T01[neighbourind]...]
    # 
    #  nc_neighbour = Array{typeof(mdbm.ncubes[1])}(undef,0)
    #  NumofNCubes=length(ncubes)
    #  for iT in 1:length(T101)
    #      Base.append!(nc_neighbour,deepcopy(ncubes))
    #      for nci in ((1+NumofNCubes*(iT-1)):(NumofNCubes+NumofNCubes*(iT-1)))
    #          nc_neighbour[nci].corner[:]=nc_neighbour[nci].corner+T101[iT].*nc_neighbour[nci].size
    #      end
    #  end
    # #-------all faceing/cornering neightbours----------------


    #-------direcational neightbours----------------
    nc_neighbour = Array{typeof(mdbm.ncubes[1])}(undef, 0)
    Nface = 2^N
    indpos = [[mdbm.T01[i][dir] for i in 1:Nface] for dir in 1:N]
    indneg = [[!mdbm.T01[i][dir] for i in 1:Nface] for dir in 1:N]
    @inbounds for nc in ncubes
        fcvals = getcornerval(nc, mdbm)
        for dir in 1:N
            # indpos=[mdbm.T01[i][dir] for i in 1:Nface]
            # indneg=[!mdbm.T01[i][dir] for i in 1:Nface]
            if issingchange(fcvals[indpos[dir]], Nf, Nc)
                push!(nc_neighbour, deepcopy(nc))
                nc_neighbour[end].corner[dir] += nc_neighbour[end].size[dir]
            end
            if issingchange(fcvals[indneg[dir]], Nf, Nc)
                push!(nc_neighbour, deepcopy(nc))
                nc_neighbour[end].corner[dir] -= nc_neighbour[end].size[dir]
            end
        end
    end
    #-------direcational neightbours----------------

    wrap_ncube_periodic!.(nc_neighbour, Ref(mdbm.axes))

    #filter!(nc -> !(any(nc.corner .< 1) || any((nc.corner + nc.size) .> [length.(mdbm.axes)...])), nc_neighbour)#remove the overhanging ncubes
    filter!(nc -> in_range(nc, mdbm.axes), nc_neighbour)
    deleteat!(nc_neighbour, is_sorted_in_sorted(nc_neighbour, ncubes))# delete the ones, whihe were in the list onriginally
    #sort!(nc_neighbour; alg=QuickSort)# it is slower the the default method
    sort!(nc_neighbour)
    unique!(nc_neighbour)
    return nc_neighbour
end
"""
    in_range(nc::NCube{IT,FT,N,Nfc}, axes::Axes{N,AT}) -> Bool

Return `true` if, for every dimension d=1…N,
  1 ≤ nc.corner[d]  and  nc.corner[d] + nc.size[d] ≤ length(axes[d]).

In other words, this returns `false` whenever any part of `nc` “sticks out”
beyond the 1…length(axes[d]) grid on dimension d.
"""
function in_range(nc::NCube{IT,FT,N,Nfc}, axes::Axes{N,AT}) where {IT,FT,N,Nfc,AT}
    @inbounds for d in 1:N
        # 1) If the cube’s “origin” is below 1 in dimension d, it’s invalid:
        if nc.corner[d] < 1
            return false
        end
        # 2) If the cube’s “far end” (corner[d] + size[d]) exceeds the number of ticks:
        if nc.corner[d] + nc.size[d] > length(axes[d])
            return false
        end
    end
    return true
end

"""
    wrap_ncube_periodic!(nc::NCube{IT,FT,N,Nfc}, axes::Axes{N,AT})

Modify `nc.corner` in‐place so that each index is wrapped into `1: length(axes[d])`.  
Equivalently, for each dimension d:
  nc.corner[d] ↦ mod(nc.corner[d] - 1, length(axes[d])) + 1

This ensures that any out‐of‐bounds corner (e.g. 0 or > length) “wraps around” periodically.
Returns `nc` for convenience.
"""
function wrap_ncube_periodic!(nc::NCube{IT,FT,N,Nfc}, axes::Axes{N,AT}) where {IT,FT,N,Nfc,AT}
    @inbounds for d in 1:N
        if axes[d].periodic
            L = length(axes[d]) - 1   # number of ticks along dimension d
            # shift corner[d] into 0..L-1 via (corner-1), then mod, then back to 1..L
            nc.corner[d] = mod(nc.corner[d] - 1, L) + 1
        end
    end
    return nc
end

"""
    checkneighbour!(mdbm::MDBM_Problem{fcT,N,Nf,Nc,t01T,t11T,IT,FT,aT}; interpolationorder::Int=0, maxiteration::Int=0)

Check the face-neighbouring n-cubes to discover the lost (missing) part of the solutions.

It is possible, that the sub-manifold (e.g.: curve) of the solution is not completely detected (it is interrupted).
This way the missing part can be explored similarly to a continuation method (with minimal function evaluation).

It works only if the dimension of the solution object is larger the zero (the maifold is not a set of poins).

    # Arguments
    - `interpolationorder::Int=0`: interpolation order method of the neighbours checked
    - `maxiteration::Int=0~: the max number of steps in the 'continuation-like' exploring. If zero, then infinity steps are allowed
"""
function checkneighbour!(mdbm::MDBM_Problem{fcT,N,Nf,Nc,t01T,t11T,IT,FT,aT}; interpolationorder::Int=0, maxiteration::Int=0, normp=20.0, ncubetolerance=0.5, doThreadprecomp=true) where {fcT,N,Nf,Nc,t01T,t11T,IT,FT,aT}#only for unite size cubes
    if isempty(mdbm.ncubes)
        println("There is no bracketing n-cubes to check!")
    else
        ncubes2check = mdbm.ncubes
        numberofiteration = 0
        while !isempty(ncubes2check) && (maxiteration == 0 ? true : numberofiteration < maxiteration)

            numberofiteration += 1
            ncubes2check = generateneighbours(ncubes2check, mdbm)
            deleteat!(ncubes2check, is_sorted_in_sorted(ncubes2check, mdbm.ncubes))#delete the ones which is already presented
            #@time fast_diff_sorted!(ncubes2check, mdbm.ncubes)#delete the ones which is already presented
            _interpolate!(ncubes2check, mdbm, Val{interpolationorder}, normp=normp, ncubetolerance=ncubetolerance, doThreadprecomp=doThreadprecomp) #remove the non-bracketing, only proper new bracketing cubes remained
            Base.append!(mdbm.ncubes, deepcopy(ncubes2check))
            sort!(mdbm.ncubes; alg=QuickSort)
        end
    end
end


"""
    connect(mdbm::MDBM_Problem{fcT,N,Nf,Nc,t01T,t11T,IT,FT,aT})

Provide the edge connection (as a list of point-paris) of the solution point-cloud based on the face-neighbouring n-cubes.
"""
function connect(mdbm::MDBM_Problem{fcT,N,Nf,Nc,t01T,t11T,IT,FT,aT}) where {fcT,N,Nf,Nc,t01T,t11T,IT,FT,aT}
    #---------- line connection (no corener neighbour is needed) --------------
    DT1 = Array{Tuple{Int64,Int64}}(undef, 0)
    for inc in 1:length(mdbm.ncubes)
        ncneigh = generateneighbours([mdbm.ncubes[inc]], mdbm)
        indinresults = index_sorted_in_sorted(ncneigh, mdbm.ncubes)#delete the ones which is already presented
        Base.append!(DT1, [(inc, x) for x in indinresults if x != 0])
    end
    filter!(d -> (d[2] > d[1]), DT1)#TODO: eleve csak ez egyik irányban levő szomszédokat kellene keresni!!!
    sort!(DT1; alg=QuickSort)
    unique!(DT1)
    return DT1
end

"""
    triangulation(DT1::Array{Tuple{Int64,Int64}})::Array{Tuple{Int64,Int64,Int64}}

Provide the triangulation (as a list of point-triplets) of the solution point-cloud based on the face-neighbouring n-cubes.
DT1 is the result of the edge connection performed by `connect`.
"""
function triangulation(DT1::Array{Tuple{Int64,Int64}})::Array{Tuple{Int64,Int64,Int64}}

    #DT=sort!()[DT1;[(d[2],d[1]) for d in DT1]]);
    DT = sort!([DT1; [d[[2, 1]] for d in DT1]])#both direction of line connection is necessay!

    #L=[filter(d->d[1]==i,DT) for i in 1:length(mdbm.ncubes)]
    L = [
        [dd[2] for dd in filter(d -> d[1] == i, DT)] for i in 1:maximum(DT1)[2]]

    DT4 = Array{Tuple{Int64,Int64,Int64,Int64}}(undef, 0)#quadratic patch
    @inbounds for i in 1:size(L, 1)
        for j in L[i]#filter((L{i}>i) for i in )
            if i > j#i must be the larges value to remove the repetition of the some surface segment
                for k in L[j]#(L{j}>i) it seems to be slower
                    if i > k# there is no backstep, and i must be the largest   (Equivalent: %((k~=i) && (i>k))%back step is not allowed
                        for m in L[k]
                            if ((m != j) && (j > m)) #&& (i>m) #back step is not allowed, and i must be the largest
                                if any(i .== L[m])
                                    push!(DT4, (i, j, k, m))
                                end
                            end
                        end
                    end
                end
            end
        end
    end

    return [[dt[[1, 2, 3]] for dt in DT4]; [dt[[3, 1, 4]] for dt in DT4]]#DT2 triangular patch from the quadratic patch

end

"""
    solve!(mdbm::MDBM_Problem{fcT,N,Nf,Nc,t01T,t11T,IT,FT,aT}, iteration::Int; interpolationorder::Int=1)

Refine the `MDBM_Problem` `iteration` times, then perform a neighbour check.
`interpolationorder` defines the interpolation method within the n-cubes.

# Examples
```julia
julia> solve!(mymdbm,4)
```
"""
function solve!(mdbm::MDBM_Problem{fcT,N,Nf,Nc,t01T,t11T,IT,FT,aT}, iteration::Int; interpolationorder::Int=1, normp=20.0, ncubetolerance=0.5, verbosity=1, doThreadprecomp=true) where {fcT,N,Nf,Nc,t01T,t11T,IT,FT,aT}
    if verbosity > 0
        @time begin
            print("  interpolate! :  ")
            @time interpolate!(mdbm, interpolationorder=interpolationorder, normp=normp, ncubetolerance=ncubetolerance, doThreadprecomp=doThreadprecomp)
            for k = 1:iteration
                print("  $k refine! & interpolate! :  ")
                @time begin
                    refine!(mdbm)
                    interpolate!(mdbm, interpolationorder=interpolationorder, normp=normp, ncubetolerance=ncubetolerance, doThreadprecomp=doThreadprecomp)
                end
            end
            print("           checkneighbour! :  ")
            @time checkneighbour!(mdbm, interpolationorder=interpolationorder, normp=normp, ncubetolerance=ncubetolerance, doThreadprecomp=doThreadprecomp)
            print("              interpolate! :  ")
            @time interpolate!(mdbm, normp=normp, ncubetolerance=ncubetolerance, doThreadprecomp=doThreadprecomp)
            print("total time:  :  ")
        end
        return mdbm
    else
        interpolate!(mdbm, interpolationorder=interpolationorder, normp=normp, ncubetolerance=ncubetolerance, doThreadprecomp=doThreadprecomp)
        for k = 1:iteration
            refine!(mdbm)
            interpolate!(mdbm, interpolationorder=interpolationorder, normp=normp, ncubetolerance=ncubetolerance, doThreadprecomp=doThreadprecomp)
        end
        checkneighbour!(mdbm, interpolationorder=interpolationorder, normp=normp, ncubetolerance=ncubetolerance, doThreadprecomp=doThreadprecomp)
        interpolate!(mdbm, normp=normp, ncubetolerance=ncubetolerance, doThreadprecomp=doThreadprecomp)
        return mdbm
    end
end

function index_sorted_in_sorted(a::AbstractVector, b::AbstractVector)::Array{Int64,1}
    # index of a[i] in b (0 if not present)
    # a and b must be sorted
    containingindex = zeros(Int64, length(a))
    startindex = 1
    @inbounds for ind2check in 1:length(a)
        detectedrange = searchsorted(b[startindex:end], a[ind2check])
        startindex = max(detectedrange.stop, detectedrange.start) + startindex - 1
        containingindex[ind2check] = isempty(detectedrange) ? 0 : startindex
        if startindex > length(b)
            break
        end
    end
    return containingindex
end

function is_sorted_in_sorted(a::AbstractVector{T}, b::AbstractVector{T})::Array{Bool,1} where T
    n = length(a)
    result = falses(n)

    i = 1           # pointer in a
    j = 1           # pointer in b
    nb = length(b)

    while i ≤ n && j ≤ nb
        ai = a[i]
        bj = b[j]

        if bj < ai
            # b[j] is too small to ever match a[i], advance j
            j += 1
        elseif bj == ai
            # found a match
            result[i] = true
            i += 1   # move on to next a[i]
        else
            # b[j] > a[i] ⇒ a[i] is not in b, record false (it's already false)
            i += 1   # advance i, keep j where it is
        end
    end

    # any remaining a[i] (i ≤ n but j>nb) must be false (already false in `result`)
    return result
end

function is_sorted_in_sorted(aVargs::AbstractVector{AT}, bfvalarg::SortedCache{AT,Tuple{RTf,RTc}})::Array{Bool,1} where {AT,RTf,RTc}
    return haskey.(Ref(bfvalarg.data), aVargs)
end

#TODO: make it faster
#function fast_diff_sorted!(A::Vector{T}, B::Vector{T}) where T
#    i = j = 1
#    result = T[]
#
#    while i <= length(A) && j <= length(B)
#        if A[i] < B[j]
#            push!(result, A[i])
#            i += 1
#        elseif A[i] > B[j]
#            j += 1
#        else  # A[i] == B[j]
#            i += 1
#            j += 1
#        end
#    end
#
#    # Add remaining elements of A
#    Base.append!(result, A[i:end])
#    return result
#end


"""
    merge_sorted(a::Vector{T}, b::Vector{T}; lt=(x,y)->x<y) where T

Merge two sorted vectors `a` and `b` into a new sorted vector.
`lt` is your “less than” comparator on elements.
"""
function merge_sorted(a::Vector{T}, b::Vector{T}; lt=(x, y) -> x < y) where T
    na, nb = length(a), length(b)
    c = Vector{T}(undef, na + nb)
    ia = ib = ic = 1

    @inbounds while ia ≤ na && ib ≤ nb
        if lt(a[ia], b[ib])
            c[ic] = a[ia]
            ia += 1
        else
            c[ic] = b[ib]
            ib += 1
        end
        ic += 1
    end

    # copy any leftovers
    @inbounds if ia ≤ na
        c[ic:end] = a[ia:na]
    else
        c[ic:end] = b[ib:nb]
    end

    return c
end

"""
    sequential_insert!(v::Vector{T}, x::T; lt = isless) where T

Insert `x` into the sorted vector `v` in‐place, preserving sorted order (using `lt(a,b)`).
Returns the updated `v`.
"""
function sequential_insert!(v::Vector{T}, x::T; lt=isless) where T
    # find the first index where `x` should go
    idx = searchsortedfirst(v, x; lt=lt)
    if idx > length(v)
        push!(v, x)
    else
        insert!(v, idx, x)
    end
    return v
end














"""
    half_merge_with_buffer!(v::Vector{T})

Assume that `v[1:mid]` and `v[mid+1:end]` are each already sorted,
where `mid = fld(length(v), 2)`.  This will merge them so that
`v[1:end]` becomes fully sorted, using a single auxiliary buffer
of size `length(v)`—no recursion, no risk of StackOverflow.

After this call, `v` is sorted.
"""
function half_merge_with_buffer!(v::Vector{T}) where T
    n = length(v)
    mid = fld(n, 2)

    # If the two runs are already disjoint (largest of first half ≤ smallest of second),
    # we’re done immediately:
    if v[mid] ≤ v[mid+1]
        return v
    end

    # 1) Allocate a temporary buffer once (same element type and length as v)
    tmp = similar(v)

    # 2) Standard merge of two sorted runs into tmp:
    i = 1          # pointer in run1  (1 ≤ i ≤ mid)
    j = mid + 1    # pointer in run2  (mid+1 ≤ j ≤ n)
    k = 1          # write‐pointer in tmp

    # Merge until one run is exhausted:
    while i ≤ mid && j ≤ n
        if v[i] ≤ v[j]
            tmp[k] = v[i]
            i += 1
        else
            tmp[k] = v[j]
            j += 1
        end
        k += 1
    end

    # Copy any leftovers from run1:
    while i ≤ mid
        tmp[k] = v[i]
        i += 1
        k += 1
    end

    # Copy any leftovers from run2:
    while j ≤ n
        tmp[k] = v[j]
        j += 1
        k += 1
    end

    # 3) Copy the merged result back into v:
    @inbounds for idx in 1:n
        v[idx] = tmp[idx]
    end

    return v
end
