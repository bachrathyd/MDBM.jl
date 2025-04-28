5+5

using GLMakie


using StaticArrays
using LinearAlgebra



N = length(mymdbm.axes)
T = Float64
pathall = Vector{SVector{N,T}}()

#for nc in mymdbm.ncubes
    k=0
    k=k+1
    nc= mymdbm.ncubes[k]


    nccorner = nc.corner
    ncsize = nc.size

    fixed_dims = Vector{Tuple{Int,Bool}}(undef, 0)
    nc_template = mymdbm.T01


    posall_tree = MDBM.PositionTree(NaN * Vector{Float64}(undef, N))
    MDBM.interpsubcubesolution!(posall_tree, [nc_template], [fixed_dims], nccorner, ncsize, mymdbm)
    #    
    #    @show posall_tree
    path2points = MDBM.extract_paths(posall_tree)
    



    #f = Figure()
    puniq=unique(collect(Iterators.flatten(path2points)))[2:end]

@show    result = linear_fit_quality(puniq)
    

#scatter!(getindex.(puniq,1),getindex.(puniq,2),getindex.(puniq,3))
scatter!(getindex.(puniq,1),getindex.(puniq,2))



    p in path2points
    #plot each line element separately
    #[plot!(collect(eachrow(reduce(hcat, pathloc)))...,lw=5) for pathloc in path2points]


    pathall_ncube = Vector{SVector{N,T}}()
    #pathloc=path2points[1]
    for pathloc in path2points
        if length(pathloc) > 2
            append!(pathall, pathloc)
            #append!(pathall_ncube, vcat(pathloc, [pathloc[2]])) #closed traingle for line plotting
        end

    end


    # #plot each ncubes separately
    # if length(pathall_ncube) > 1
    #     plot!(collect(eachrow(reduce(hcat, pathall_ncube)))[1:3]..., lw=5)
    # end
    append!(pathall, pathall_ncube)
#end

faces=(1:4:(length(pathall))) .+ collect(1:N)'
mesh!(pathall,faces,alpha=1.0)

