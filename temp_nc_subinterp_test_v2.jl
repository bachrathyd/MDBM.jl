5 + 5

using GLMakie


using StaticArrays
using LinearAlgebra

begin

N = length(mymdbm.axes)
T = Float64
pathall = Vector{SVector{N,T}}()
error_all = Vector{Float64}()

for nc in mymdbm.ncubes
    # k = 0
    # k = k + 1
    # nc = mymdbm.ncubes[k]


    nccorner = nc.corner
    ncsize = nc.size

    fixed_dims = Vector{Tuple{Int,Bool}}(undef, 0)
    nc_template = mymdbm.T01

    
    posall_tree = MDBM.PositionTree(nc.posinterp)
    edges, all_edge_fixed_dims = MDBM.generate_sub_faces(nc_template, fixed_dims)

    MDBM.interpsubcubesolution!(posall_tree, edges, all_edge_fixed_dims, nccorner, ncsize, mymdbm)

    path2points = MDBM.extract_paths(posall_tree)




    puniq = unique(collect(Iterators.flatten(path2points)))

    puniq_pos=getinterpolatedsolution.(puniq, Ref(nc.corner), Ref(mymdbm.axes))
    path2points_pos = [getinterpolatedsolution.(inner, Ref(nc.corner), Ref(mymdbm.axes)) for inner in path2points]

    xyz_sol = nc.posinterp
    gxyz = nc.gradient

    total = 0.0
    @inbounds for g in eachcol(gxyz)
        for p in puniq
            total += abs(dot(p .- xyz_sol, g) / norm(g))
        end
    end

    #    arrowcolor = strength, linecolor = strength)

  #  if N == 3
  #      scatter!(getindex.(puniq_pos, 1), getindex.(puniq_pos, 2), getindex.(puniq_pos, 3))
  #  else
  #      scatter!(getindex.(puniq_pos, 1), getindex.(puniq_pos, 2))
  #  end


    #plot each line element separately
    #[plot!(collect(eachrow(reduce(hcat, pathloc)))...,lw=5) for pathloc in path2points]

    # pathall = Vector{SVector{N,T}}()


    pathall_ncube = Vector{SVector{N,T}}()
    #pathloc=path2points[1]
    for pathloc in path2points_pos
        if length(pathloc) > 2

            push!(pathall, fill(NaN, N))
            append!(pathall, pathloc)
            push!(error_all, NaN)
            append!(error_all, ones(length(pathloc)) .* total)
            #append!(pathall_ncube, vcat(pathloc, [pathloc[2]])) #closed traingle for line plotting
        end
    end
    #lines!(getindex.(pathall,1),getindex.(pathall,2),color=getindex.(pathall,1) .* 0 .+ total)
    #lines!(getindex.(pathall,1),getindex.(pathall,2),getindex.(pathall,3),color=getindex.(pathall,1) .* 0 .+ total)

    # #plot each ncubes separately
    # if length(pathall_ncube) > 1
    #     plot!(collect(eachrow(reduce(hcat, pathall_ncube)))[1:3]..., lw=5)
    # end
    #append!(pathall, pathall_ncube)
end

end


f = Figure()
n_faces = (1:4:(length(pathall))) .+ collect(1:N)'
mesh(pathall, n_faces, alpha=1.0)


if N == 3
    lines!(getindex.(pathall, 1), getindex.(pathall, 2), getindex.(pathall, 3), color=log.(error_all), linewidth=3)
else
    lines!(getindex.(pathall, 1), getindex.(pathall, 2), color=log.(error_all), linewidth=3)
end
#plot(getindex.(pathall,1),getindex.(pathall,2),getindex.(pathall,3),lw=5)


Colorbar(f[1, 2], limits=(log(minimum(error_all)), log(maximum(error_all))), colormap=:viridis,    flipaxis=false)

f

lines(getindex.(pathall, 1), getindex.(pathall, 2), getindex.(pathall, 3), color=log.(error_all), linewidth=3)

