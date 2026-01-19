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
# function foo_par2_codim1(x, y)
#     pow = 3
#     abs(x)^pow + abs(y)^pow - 2.0^pow
# end
# 
# mymdbm = MDBM_Problem(foo_par2_codim1, [-3.0:0.5:3.0, -3.0:0.5:3.0])
# solve!(mymdbm, 2, interpolationorder=1, ncubetolerance=0.0)

function foo_par2_codim1(x, y)
    norm([x,y],110.9999)-1.0#((x^2.0 + y) - 1.0^2.0)*x #TODO: test with this function, too
    #((x^2.0 + y) - 1.0^2.0)*x #TODO: test with this function, too
end

# refines only n-cubes where the error is greater than 50% betweenthe worst and best error
mymdbm = MDBM_Problem(foo_par2_codim1, [-3.1:3.0, -3.1:3.0])
@time solve!(mymdbm, 1,refinementratio=0.5)#7)#number of refinements - increase it
## Itrative refinement ------------------------------------

for _ in 1:2
    nc_list = 1:size(mymdbm.ncubes, 1)


    errv_s = [MDBM.getscaled_local_point(MDBM.ncube_error_vector(nc), nc, mymdbm.axes) for nc in mymdbm.ncubes]
    err_norm = norm.(errv_s)
    # nc_list = nc_list[err_norm.>=sort(err_norm)[length(err_norm)÷10*7+1]]
    nc_list = nc_list[err_norm.>=(minimum(err_norm)+maximum(err_norm))/2]

    #    nc_list= nc_list[err_norm .> 0.0001]


    # TODO: this is a problem, doubling must be done for only he cubes which hase size 1, and at the location where it is size 1 - this way i will not be able to tell the size (diffference) of the neighbouring n-cubes
    nc_size_minimum_in_the_list = minimum([minimum(nc.size) for nc in mymdbm.ncubes[nc_list]])
    if nc_size_minimum_in_the_list < 2
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
    MDBM.interpolate!(mymdbm, interpolationorder=1, normp=10.0, ncubetolerance=ncubetolerance = 0.4)
    #end
    #@time solve!(mymdbm, 1,interpolationorder=1)


         #        ncube_layers[]  = [    
         #            mymdbm.ncubes,
         #   ]

end


##-----------------------





# 2. Set up the interactive plot
GLMakie.closeall()
GLMakie.activate!(; title="Interactive Neighbour Check")
fig = Figure(size=(1200, 1000))
ax = GLMakie.Axis(fig[1, 1], title="Click on a cube to see its neighbours", aspect=DataAspect())

# 3. Define the layer-based plotting data structure and colors
# The central data structure is a list of layers, where each layer is a list of n-cubes.
ncube_layers = Observable{Vector{Vector{MDBM.NCube}}}([mymdbm.ncubes])

# Define a color for each layer. The format is (color, alpha).
layer_colors = [
    (:lightblue, 0.5), # Layer 1: All cubes
    (:red, 0.5),       # Layer 2: Selected cube
    (:green, 0.5),     # Layer 3: Neighbours
    (:purple, 0.5),    # Layer 4: (Placeholder for filtered neighbours)
    (:orange, 0.5)     # Layer 5... and so on
]

# 4. Set up the plotting logic to redraw everything when ncube_layers changes
on(ncube_layers) do layers
    empty!(ax) # Clear the axis before redrawing
    for (i, ncube_list) in enumerate(layers)
        # Skip empty layers
        if isempty(ncube_list)
            continue
        end

        # Get color for the current layer, cycling through if needed
        color = layer_colors[mod1(i, length(layer_colors))]

        # Create polygons and plot them
        polys = [Polygon(get_cube_corners(nc, mymdbm)) for nc in ncube_list]
        poly!(ax, polys, color=color, strokecolor=:black, strokewidth=1)
    end
end

# Trigger initial plot
notify(ncube_layers)


# 5. Create a lookup map for faster/more reliable neighbour index finding
ncube_to_idx = Dict(nc => i for (i, nc) in enumerate(mymdbm.ncubes))

# 6. Register the mouse click event
on(events(ax).mousebutton, priority=2) do event
    
@time solve!(mymdbm, 1,refinementratio=0.5)#7)#number of refinements - increase it
    if event.button == Mouse.left && event.action == Mouse.press
        # Get mouse position in data coordinates
        mx, my = mouseposition(ax)

        # Find the clicked n-cube
        clicked_idx = find_clicked_ncube_idx(Point2f(mx, my), mymdbm)


        if clicked_idx !== nothing
            println("Selected n-cube index: ", clicked_idx)

            selected_ncube = mymdbm.ncubes[clicked_idx]
            neighbours = MDBM.generateneighbours([selected_ncube], mymdbm)
            neighbours = MDBM.generateneighbours(mymdbm.ncubes, mymdbm)

            neighbours_filter = deepcopy(neighbours)
            #  filter!(x -> !(x in mymdbm.ncubes), neighbours_filter)
            filter!((nc) -> !MDBM.is_overlapping(nc, mymdbm.ncubes), neighbours_filter)
            # Define the layers to plot
            # Layer 1: All cubes
            # Layer 2: The selected cube
            # Layer 3: The neighbours
         
            new_layers = [
                mymdbm.ncubes,
                [selected_ncube],
                #neighbours,
                neighbours_filter,
            ]
            # 
            # Update the observable to trigger a replot
            ncube_layers[] = new_layers

        else
            println("No cube selected. Clearing selection.")
            # Reset to only showing the base layer
            ncube_layers[] = [mymdbm.ncubes]
        end
    end
    # Let other events pass through
    return Consume(false)
end

# Display the figure
display(fig)

println("\nSetup complete. Click on a blue cube in the plot to highlight it (red) and its neighbours (green).")
# 
# 
# ##
# solve!(mymdbm, 2, interpolationorder=1, ncubetolerance=0.0)
# 
#                  ncube_layers[]  = [    
#                      mymdbm.ncubes,
#             ]
# ##
# mymdbm_start=deepcopy(mymdbm);
# ##
# mymdbm=deepcopy(mymdbm_start);
# #
# mdbm=deepcopy(mymdbm)
# mymdbm=deepcopy(mdbm)
# empty!(mymdbm.ncubes)
# append!(mymdbm.ncubes,deepcopy(mdbm.ncubes[1:30]))
# 
#     ncube_layers[]  = [    
#                      mymdbm.ncubes,
#             ];
# 
# @time checkneighbour!(mymdbm,verbosity=0);
# @profview  checkneighbour!(mymdbm,verbosity=0);
#     ncube_layers[]  = [    
#                      mymdbm.ncubes,
#             ];
# 
# 
# ##
# 
# 