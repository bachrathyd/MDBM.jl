using Revise
using MDBM
#using GLMakie
#GLMakie.closeall()
#GLMakie.activate!(; title="3 parameters, codimension 1")
#-----------------------------

#This type of type definitsion is an overkill here!
function foo_par3_codim1(x::Float64, y::Float64, z::Float64)::Float64
    x^2.0 + y^2.0 + z^2.0 - 2.0^2.0::Float64
end

# constraint - calculate only the points where the constraint is satisfied (e.g.: on the positiev side)
function c(x::Float64, y::Float64, z::Float64)::Float64
    x^2.0 + y^2.0 - 0.5^2.0
end


mymdbm = MDBM_Problem(foo_par3_codim1, [-3.0:1.0, -1.0:3.0, -1.0:3.0], constraint=c)
solve!(mymdbm, 2, doThreadprecomp=true, verbosity=0)




function foo(n, b=true)
    mymdbm = MDBM_Problem(foo_par3_codim1, [-3.0:1.0, -1.0:3.0, -1.0:3.0], constraint=c, memoization=b)
    solve!(mymdbm, n, doThreadprecomp=true, verbosity=0)
end

for kk in 1:6
    println(kk)
    @time foo(kk, false)#number of refinements - increase it slightly to see smoother results 
end

##
5+5
#foo(2,false)
#@profview foo(4);

##  using BenchmarkTools
##  @benchmark foo($5)
##  f = Figure(size=(1000, 600))
##  ax1 = GLMakie.Axis3(f[1, 1])
##  
##  xyz_sol = getinterpolatedsolution(mymdbm)
##  scatter!(ax1, xyz_sol..., markersize=15, color=:red, marker='x', strokewidth=3, label="solution")
##  
##  # show the points where the function is evaluated
##  xyz_val = getevaluatedpoints(mymdbm)
##  fval = getevaluatedfunctionvalues(mymdbm)
##  scatter!(xyz_val..., color=sign.(fval), label="evaluated")
##  
##  # connecting and plotting the "mindpoints" of the n-cubes
##  DT1 = connect(mymdbm)
##  edge2plot_xyz = [reduce(hcat, [i_sol[getindex.(DT1, 1)], i_sol[getindex.(DT1, 2)], fill(NaN, length(DT1))])'[:] for i_sol in xyz_sol]
##  lines!(edge2plot_xyz..., linewidth=5, label="midpoints solution connected")
##  
##  # #plotting the gradintes
##  gxyz = getinterpolatedgradient(mymdbm.ncubes, mymdbm)
##  arrows!(xyz_sol..., gxyz[1]..., arrowsize=0.1, lengthscale=0.1, label="gradient")#    arrowcolor = strength, linecolor = strength)
##  
##  
##  
##  
##  
##  #--------------------------- Sub-cube interpolation----------------
##  ax2 = GLMakie.Axis3(f[1, 2])
##  #calcuatin the sub-cubes interpolations stored in the mymdbm.ncubes[i].posinterp
##  interpsubcubesolution!(mymdbm)
##  #extracting the resutls to from the 
##  path2points = extract_paths(mymdbm)
##  
##  #extracting the unique points and plotting
##  puniq = unique(collect(Iterators.flatten(Iterators.flatten(path2points))))
##  scatter!(ax2, getindex.(puniq, 1), getindex.(puniq, 2), getindex.(puniq, 3), markersize=10, color=:green, label="subface - solution")
##  
##  
##  
##  #exctracing the simplexes for each ncube
##  flatened_path2points = collect(Iterators.flatten(path2points))
##  #eliminating the points with less than 2 points (caused by fininte precision)
##  true_truflatened_path2points = flatened_path2points[length.(flatened_path2points).==3]
##  
##  
##  #plotting the lines between the points
##  n_faces = reshape(1:(3*length(true_truflatened_path2points)), (3, length(true_truflatened_path2points)))'
##  vertices_mat = hcat(Iterators.flatten(true_truflatened_path2points)...)
##  mesh!(vertices_mat, n_faces, alpha=0.5, label="subface - local simplex")
##  
##  
##  display(GLMakie.Screen(), f)
##  axislegend(ax1)
##  axislegend(ax2)
##  
##  
##  
##  #Final plots with diffent methods
##  ff = Figure(size=(1000, 600))
##  lines(ff[1, 1], edge2plot_xyz..., linewidth=5)
##  mesh(ff[1, 2], vertices_mat, n_faces; color=:steelblue, shininess=3000)
##  supertitle = Label(ff[0, :], "Final plots", fontsize=20)
##  display(GLMakie.Screen(), ff)
##  #Fun with lights
##  
##  fig = Figure(size=(1000, 900))
##  ax11 = LScene(fig[1, 1], scenekw=(lights=[DirectionalLight(RGBf(0, 0, 0), Vec3f(-1, 0, 0))],))
##  ax12 = LScene(fig[1, 2], scenekw=(lights=[DirectionalLight(RGBf(1, 1, 1), Vec3f(-1, 0, 0))],))
##  lights = [
##      DirectionalLight(RGBf(0, 0, 0.7), Vec3f(-1, -1, 0)),
##      DirectionalLight(RGBf(0.7, 0.2, 0), Vec3f(-1, 1, -1)),
##      DirectionalLight(RGBf(0.7, 0.7, 0.7), Vec3f(1, -1, -1))
##  ]
##  ax21 = LScene(fig[2, 1], scenekw=(lights=lights,))
##  ax22 = LScene(fig[2, 2], scenekw=(lights=[DirectionalLight(RGBf(4, 2, 1), Vec3f(0, 0, -1))],))
##  for ax in (ax11, ax12, ax21, ax22)
##      mesh!(ax, vertices_mat, n_faces; color=:steelblue)
##  end
##  supertitle = Label(fig[0, :], "Fun with lights", fontsize=20)
##  display(GLMakie.Screen(), fig)
##  
##  
##  
##  ##
##  
##  
##  using LinearAlgebra
##  function export_mesh_to_stl(filename::String, vertices_mat::AbstractMatrix, n_faces::AbstractMatrix)
##      open(filename, "w") do io
##          println(io, "solid mesh")
##  
##          for j in 1:size(n_faces, 2)
##              # Get indices of triangle vertices (assuming 0-based indices from GLMakie)
##              i1, i2, i3 = n_faces[:, j]
##  
##              # Get triangle vertices
##              v1 = vertices_mat[:, i1]
##              v2 = vertices_mat[:, i2]
##              v3 = vertices_mat[:, i3]
##  
##              # Compute normal (normalized cross product)
##              n = cross(v2 - v1, v3 - v1)
##              norm_n = norm(n)
##              n = norm_n == 0 ? [0.0, 0.0, 0.0] : n / norm_n
##  
##              # Write triangle
##              println(io, "  facet normal $(n[1]) $(n[2]) $(n[3])")
##              println(io, "    outer loop")
##              println(io, "      vertex $(v1[1]) $(v1[2]) $(v1[3])")
##              println(io, "      vertex $(v2[1]) $(v2[2]) $(v2[3])")
##              println(io, "      vertex $(v3[1]) $(v3[2]) $(v3[3])")
##              println(io, "    endloop")
##              println(io, "  endfacet")
##          end
##  
##          println(io, "endsolid mesh")
##      end
##  end
##  
##  
##  using LinearAlgebra
##  
##  function export_mesh_to_binary_stl(filename::String, vertices::AbstractMatrix, faces::AbstractMatrix)
##      n_tri = size(faces, 2)
##      open(filename, "w") do io
##          # 80-byte header
##          hdr = zeros(UInt8, 80)
##          txt = "Binary STL from Julia"
##          hdr[1:length(txt)] = codeunits(txt)
##          write(io, hdr)
##  
##          # write triangle count (UInt32 little-endian)
##          write(io, UInt32(n_tri))
##  
##          for j in 1:n_tri
##              i1, i2, i3 = faces[:, j]
##              v1, v2, v3 = vertices[:, i1], vertices[:, i2], vertices[:, i3]
##              n = cross(v2 - v1, v3 - v1)
##              if norm(n) > 0
##                  n ./= norm(n)
##              else
##                  n .= 0.0
##              end
##  
##              # write normal + verts as Float32 LE
##              write(io,
##                  Float32(n[1]), Float32(n[2]), Float32(n[3]),
##                  Float32(v1[1]), Float32(v1[2]), Float32(v1[3]),
##                  Float32(v2[1]), Float32(v2[2]), Float32(v2[3]),
##                  Float32(v3[1]), Float32(v3[2]), Float32(v3[3]),
##                  UInt16(0)           # attribute byte count
##              )
##          end
##      end
##  end
##  
##  export_mesh_to_stl("output_h_sqrt2_tau.stl", vertices_mat, n_faces')
##  export_mesh_to_binary_stl("output_binary_h_sqrt2_tau.stl", vertices_mat, n_faces')