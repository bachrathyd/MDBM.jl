using Revise
using MDBM
using GLMakie
using LinearAlgebra
using GeometryBasics

# Helper function to get the corner coordinates of an n-cube
function get_cube_corners(nc, mdbm)
    corners_for_nc = MDBM.corner([nc], mdbm.T01)[1]
    
    # The user noted the correct order to draw a 2D square is [1, 2, 4, 3]
    correct_order_indices = [1, 2, 4, 3]
    ordered_corners = corners_for_nc[correct_order_indices]
    
    # Convert corner indices to data coordinates
    points = [Point2f(mdbm.axes[1][c[1]], mdbm.axes[2][c[2]]) for c in ordered_corners]
    
    return points
end

# Helper function to find which n-cube contains the clicked point
function find_clicked_ncube_idx(pos, mdbm)
    for (i, nc) in enumerate(mdbm.ncubes)
        # Get the min/max corners of the n-cube in data coordinates
        min_corner_idx = nc.corner
        max_corner_idx = nc.corner .+ nc.size
        min_coord = [mdbm.axes[d][min_corner_idx[d]] for d in 1:length(pos)]
        max_coord = [mdbm.axes[d][max_corner_idx[d]] for d in 1:length(pos)]
        
        # Check if the point is inside the bounding box
        if all(min_coord .<= pos .<= max_coord)
            return i
        end
    end
    return nothing
end

# 1. Define the problem and solve it to get a set of n-cubes
function foo_par2_codim1(x, y)
    pow = 3
    abs(x)^pow + abs(y)^pow - 2.0^pow
end

mymdbm = MDBM_Problem(foo_par2_codim1, [-3.0:0.5:3.0, -3.0:0.5:3.0])
solve!(mymdbm, 0, interpolationorder=1, ncubetolerance=0.0)



for _ in 1:3  
     nc_list = 1:size(mymdbm.ncubes, 1)


    errv_s = [MDBM.getscaled_local_point(MDBM.ncube_error_vector(nc), nc, mymdbm.axes) for nc in mymdbm.ncubes]
    err_norm = norm.(errv_s)
   # nc_list = nc_list[err_norm.>=sort(err_norm)[length(err_norm)÷10*7+1]]
    nc_list = nc_list[err_norm.>=(minimum(err_norm)+maximum(err_norm))/2]

   #    nc_list= nc_list[err_norm .> 0.0001]


    # TODO: this is a problem, doubling must be done for only he cubes which hase size 1, and at the location where it is size 1 - this way i will not be able to tell the size (diffference) of the neighbouring n-cubes
     nc_size_minimum_in_the_list = minimum([minimum(nc.size) for nc in mymdbm.ncubes[nc_list]])
     if nc_size_minimum_in_the_list<2
            MDBM.doubling!(mymdbm, [1, 2])
     end
     
    if length(nc_list) == 0
        break
    end
    @show length(nc_list)
    #MDBM.refinencubes!(mymdbm.ncubes, 1:size(mymdbm.ncubes, 1) ÷ 2, [1, 2])
    MDBM.refinencubes!(mymdbm.ncubes, nc_list, [1, 2])
    #MDBM.refinencubes!(mymdbm.ncubes, 1:size(mymdbm.ncubes, 1) , [1])
    #MDBM.refinencubes!(mymdbm.ncubes, 1:size(mymdbm.ncubes, 1) , [2])
    MDBM.interpolate!(mymdbm, interpolationorder=1,normp=10.0, ncubetolerance=ncubetolerance=0.4)
    #end
    #@time solve!(mymdbm, 1,interpolationorder=1)



end
##





# 2. Set up the interactive plot
GLMakie.closeall()
GLMakie.activate!(; title="Interactive Neighbour Check")
fig = Figure(size=(1200, 1000))
ax = GLMakie.Axis(fig[1, 1], title="Click on a cube to see its neighbours", aspect=DataAspect())

# Observables to store the state
selected_idx = Observable{Union{Int, Nothing}}(nothing)
neighbour_indices = Observable{Vector{Int}}([])

# 3. Plot all the n-cubes with a default color
all_ncubes_polys = [Polygon(get_cube_corners(nc, mymdbm)) for nc in mymdbm.ncubes]
poly!(ax, all_ncubes_polys, color=:lightblue, strokecolor=:black, strokewidth=1)

# 4. Define the plotting for selected and neighbour cubes
# We create a function that returns the polygons for the selected and neighbour cubes
function get_highlight_polys(selected_idx, neighbour_indices)
    selected_poly = Vector{Polygon}()
    neighbour_polys = Vector{Polygon}()

    if selected_idx !== nothing
        push!(selected_poly, Polygon(get_cube_corners(mymdbm.ncubes[selected_idx], mymdbm)))
    end
    
    for idx in neighbour_indices
        push!(neighbour_polys, Polygon(get_cube_corners(mymdbm.ncubes[idx], mymdbm)))
    end
    
    return selected_poly, neighbour_polys
end

# Use lift to create dependent observables for the polygons.
# We explicitly type the generated vectors as Vector{Polygon} so that when the
# underlying observables are empty, we get an empty Vector{Polygon} instead of
# an empty Vector{Any}, which Makie's poly() cannot handle.
selected_poly = @lift($selected_idx === nothing ? Vector{Polygon}() : Polygon[Polygon(get_cube_corners(mymdbm.ncubes[$selected_idx], mymdbm))])
neighbour_polys = @lift(Polygon[Polygon(get_cube_corners(mymdbm.ncubes[i], mymdbm)) for i in $neighbour_indices])


# Plot the highlighted polygons
poly!(ax, selected_poly, color=:red, strokecolor=:black, strokewidth=2)
poly!(ax, neighbour_polys, color=:lightgreen, strokecolor=:black, strokewidth=1)

# 5. Create a lookup map for faster/more reliable neighbour index finding
ncube_to_idx = Dict(nc => i for (i, nc) in enumerate(mymdbm.ncubes))

# 6. Register the mouse click event
on(events(ax).mousebutton, priority = 2) do event
    if event.button == Mouse.left && event.action == Mouse.press
        # Get mouse position in data coordinates
        mx, my = mouseposition(ax)
        
        # Find the clicked n-cube
        clicked_idx = find_clicked_ncube_idx(Point2f(mx, my), mymdbm)
        
        if clicked_idx !== nothing
            println("Selected n-cube index: ", clicked_idx)
            
            # Update the selected index observable
            selected_idx[] = clicked_idx
            
            # Find neighbours
            selected_ncube = mymdbm.ncubes[clicked_idx]
            neighbours = MDBM.generateneighbours([selected_ncube], mymdbm)
            
            # Get indices from the lookup map
            neighbour_idxs = [get(ncube_to_idx, nb, nothing) for nb in neighbours]
            filter!(!isnothing, neighbour_idxs)
            
            println("Found ", length(neighbour_idxs), " neighbours with indices: ", neighbour_idxs)
            
            # Update the neighbours observable
            selected_idx[] = nothing
            neighbour_indices[] = convert(Vector{Int}, neighbour_idxs)
        else
            println("No cube selected. Clearing selection.")
            selected_idx[] = nothing
            neighbour_indices[] = []
        end
    end
    # Let other events pass through
    return Consume(false)
end

# Display the figure
display(fig)

println("\nSetup complete. Click on a blue cube in the plot to highlight it (red) and its neighbours (green).")

