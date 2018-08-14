type mdbm_object
    f::Function## evaluated function
    fconstrain::Function## evaluated function
    fvectorized::Bool ## is function f can be called in a vectorized form
    ax::Array{Array{Float64,1},1}##description of the initial grid
    Ndim::Int64
    Nax::Array{Int64,1}
    Naxstride::Array{Int64,1}
  
    Ncodim::Int64
    HC::Array{Float64,2}#T computed function values and constrain value
    linindex::Array{Int64,1} # corresponding linear-indices
    vectindex::Array{Int64,2}# corresponding sub-indices
    #N::Int64
    #compind
    #pointerp
    #DT
    ncubelin::Array{Int64,1} # linear-indices of the bracketing n-cubes
    ncubevect::Array{Int64,2} # sub-indices of the bracketing n-cubes
  
    posinterp::Array{Float64,2}# the interpolated valuse within the bracketing n-cubes
    gradient::Array{Float64,3}# the corresponding gradients within the bracketing n-cubes
  
    DT1::Array{Int64,2}#line 'tiangulation' of the resultant interpolated values (basd on the n-cube sub-indices)
    DT2::Array{Int64,2}#surface tiangulation of the resultant interpolated values (basd on the n-cube sub-indices)
  
    isconstrained::Bool #is f provides contrain? all the constraints are combinded!!! length C===1
    interporder::Int64 # 0,1,2
    selectionmode::Int64 # 0-safe selection, 1st order interpolation mased,2???
  end