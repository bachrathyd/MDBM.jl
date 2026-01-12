#Attial Kossa  -- Ogedn stability example - using MDBM to find the stability boundary

5 + 5
using Revise
using MDBM
using GLMakie

using LinearAlgebra

using ForwardDiff
GLMakie.closeall()
GLMakie.activate!(; title="2 parameters, codimension 1")
#-----------------------------

# function R_Ogden(λ, λT, ν0, α)
#     term1 = (1 / α) * (λ^(-α / 3 - 1)) / (λT^(2 * α / 3 + 2)) * (λT^α - λ^α)
#     term2 = ((1 + ν0) / (1 - 2 * ν0)) * (λ * λT^2 - 1)
#     return (2 / 3) * (term1 + term2)
# end

function R_Ogden(λ, λT, ν0, α)
    term1 = (1 - 2 * ν0)*(1 / α) * (λ^(-α / 3 - 1)) / (λT^(2 * α / 3 + 2)) * (λT^α - λ^α)
    term2 = ((1 + ν0)) * (λ * λT^2 - 1)
    return (2 / 3) * (term1 + term2)
end

## Reformulated
#function R_Ogden(λ, λT, ν0, α)
#    
#        term1 =0
#        term2 = 0
#   if α => 0.0
#        term1 = (λ^(-α / 3 - 1)) * (λT^α - λ^α) * (1 - 2 * ν0)
#        term2 = α*(1 + ν0) * (λ * λT^2 - 1) * (λT^(2 * α / 3 + 2))
#   else
#        term1 = -(λ^(-α / 3 - 1)) * (λT^α - λ^α) * (1 - 2 * ν0)
#        term2 = -α*(1 + ν0) * (λ * λT^2 - 1) * (λT^(2 * α / 3 + 2))
#    end
#    return (2 / 3) * (term1 + term2)
#end

R_Ogden(0.1847, 0.0, 0.45, 0.001)
R_Ogden(0.1847, 0.0, 0.45, -0.001)

"""
Calculates the derivative of R_Ogden with respect to λT.

By default, it uses automatic differentiation (ForwardDiff).
If `dofinitediff=true` is passed, it uses a central finite difference
method with step size `h`.
"""
function dR_Ogden_dλT(λ::Float64,
    λT::Float64,
    ν0::Float64,
    α::Float64;  # Semicolon separates positional from keyword args
    dofinitediff::Bool=false,
    h::Float64=1e-2)::Float64

    if dofinitediff
        # --- Central Finite Difference Method ---
        # f'(x) ≈ (f(x+h) - f(x-h)) / (2h)        
        f_plus = R_Ogden(λ, λT + h, ν0, α)
        f_minus = R_Ogden(λ, λT - h, ν0, α)

        return (f_plus - f_minus) / (2.0 * h)

    else
        # --- Automatic Differentiation (ForwardDiff) Method ---
        # Original code
        # We create an anonymous function that only takes the parameter of interest
        return ForwardDiff.derivative(lamdaT -> R_Ogden(λ, lamdaT, ν0, α), λT)

    end
end




"""
Calculates the derivative of R_Ogden with respect to α.

By default, it uses automatic differentiation (ForwardDiff).
If `dofinitediff=true` is passed, it uses a central finite difference
method with step size `h`.
"""
function dR_Ogden_dα(λ::Float64,
    λT::Float64,
    ν0::Float64,
    α::Float64;  # Semicolon separates positional from keyword args
    dofinitediff::Bool=false,
    h::Float64=1e-2)::Float64

    if dofinitediff
        # --- Central Finite Difference Method ---
        # f'(x) ≈ (f(x+h) - f(x-h)) / (2h)

        f_plus = R_Ogden(λ, λT, ν0, α + h)
        f_minus = R_Ogden(λ, λT, ν0, α - h)

        return (f_plus - f_minus) / (2.0 * h)

    else
        # --- Automatic Differentiation (ForwardDiff) Method ---
        # Original code
        # We create an anonymous function that only takes the parameter of interest
        f = α_par -> R_Ogden(λ, λT, ν0, α_par)

        return ForwardDiff.derivative(f, α)
    end
end


const α_fix = 1.0
const v0_fix = 0.45


axis_λ = MDBM.Axis(collect(LinRange(0.05, 1.0, 5)), "λ")
axis_λT = MDBM.Axis(collect(LinRange(0.02, 4.0, 5)), "λT")
#axis_v0 = MDBM.Axis(collect(LinRange(0.1, 0.649, 5)), "v0")
axis_v0 = MDBM.Axis(collect(LinRange(0.1, 0.5, 5)), "v0")
axis_alfa = MDBM.Axis(collect(LinRange(-2.0, 5.0, 5)), "α")

dR_Ogden_dλT(0.1847, 1.4, 0.45, -1.0)


R_Ogden_λλT(λ::Float64, λT::Float64) = R_Ogden(λ, λT, v0_fix, α_fix)::Float64
dR_Ogden_dlambda_λλT(λ::Float64, λT::Float64) = dR_Ogden_dλT(λ, λT, v0_fix, α_fix)::Float64
#dR_Ogden_dlambda_λλT(λ::Float64, λT::Float64)=dR_Ogden_dλT(λ, λT, v0_fix,α_fix,dofinitediff=true)::Float64
@code_warntype R_Ogden(0.1847, 1.4, 0.45, -1.0)
@code_warntype dR_Ogden_dλT(0.1847, 1.4, 0.45, -1.0)
@code_warntype dR_Ogden_dα(0.1847, 1.4, 0.45, -1.0)


dR_Ogden_dlambda_λλT(0.1847, 1.45)
dR_Ogden_dlambda_λλT(0.1847, 1.5)
dR_Ogden_dlambda_λλT(0.1847, 1.40)


Ogden_mdbm = MDBM_Problem(R_Ogden_λλT, [axis_λ, axis_λT])
@time solve!(Ogden_mdbm, 6,interpolationorder=1)#number of refinements - increase it slightly to see smoother results 
Ogden_deiv_mdbm = MDBM_Problem(dR_Ogden_dlambda_λλT, [axis_λ, axis_λT])
@time solve!(Ogden_deiv_mdbm, 6,interpolationorder=1)#number of refinements - increase it slightly to see smoother results 



f = Figure(size=(1400, 1000))
#show the final resolution of the grid based on the minorticks
kwargs = (; xminorticksvisible=true, xminorgridvisible=true, yminorticksvisible=true, yminorgridvisible=true)
ax1 = GLMakie.Axis(f[1, 1]; xminorticks=Ogden_mdbm.axes[1].ticks, yminorticks=Ogden_mdbm.axes[2].ticks, kwargs..., xlabel="λ", ylabel="λT")

# # n-cube interpolation
xy_sol = getinterpolatedsolution(Ogden_mdbm)
# scatter!(xy_sol..., markersize = 15, color = :red,marker ='x',strokewidth=3,label = "solution")
# 
# # show the points where the function is evaluated
# xy_val = getevaluatedpoints(Ogden_mdbm)
# fval=getevaluatedfunctionvalues(Ogden_mdbm)
# scatter!(xy_val...,color=sign.(fval),label = "evaluated")

# connecting and plotting the "mindpoints" of the n-cubes
DT1 = connect(Ogden_mdbm)
edge2plot_xyz = [reduce(hcat, [i_sol[getindex.(DT1, 1)], i_sol[getindex.(DT1, 2)], fill(NaN, length(DT1))])'[:] for i_sol in xy_sol]
lines!(edge2plot_xyz..., linewidth=5, label="midpoints solution connected")

# connecting and plotting the "mindpoints" of the n-cubes
xy_sol = getinterpolatedsolution(Ogden_deiv_mdbm)
DT1 = connect(Ogden_deiv_mdbm)
edge2plot_xyz = [reduce(hcat, [i_sol[getindex.(DT1, 1)], i_sol[getindex.(DT1, 2)], fill(NaN, length(DT1))])'[:] for i_sol in xy_sol]
lines!(edge2plot_xyz..., linewidth=5, label="midpoints solution connected")


## --------------------- 3D case - 1 equtaion ------------------------


R_Ogden_λλTv0(λ::Float64, λT::Float64, v0::Float64) = R_Ogden(λ, λT, v0, α_fix)::Float64

Ogden_3D_mdbm = MDBM_Problem(R_Ogden_λλTv0, [axis_λ, axis_λT, axis_v0])
@time solve!(Ogden_3D_mdbm, 6, interpolationorder=1, checkneighbourNum=0,ncubetolerance=-0.2)#number of refinements - increase it slightly to see smoother results 

@time checkneighbour!(Ogden_3D_mdbm; interpolationorder=1, maxiteration=Inf, normp=20.0, ncubetolerance=0.25, doThreadprecomp=true)


f = Figure(size=(1400, 1000))
ax2 = GLMakie.Axis3(f[1, 1], xlabel="λ", ylabel="λT", zlabel="v0")

x, y, z = getinterpolatedsolution(Ogden_3D_mdbm)
scatter!(ax2, x, y, z, markersize=5, label="solution", color=:black)



# # one observable holding all 3D points
# pts = Observable(Point3f[])
# scatter!(ax2, pts; markersize = 3)
# 
# record(f, "ogden_scan.gif", 1:11; framerate = 5) do k
#     println("Frame $k")
#    @time  checkneighbour!(Ogden_3D_mdbm;
#         interpolationorder = 1,
#         maxiteration       = 1,
#         normp              = 20.0,
#         ncubetolerance     = 0.25,
#         doThreadprecomp    = true)
# 
#     x, y, z = getinterpolatedsolution(Ogden_3D_mdbm)
# 
#     # ensure equal length; trim to the shortest if needed
#     n = min(length(x), length(y), length(z))
#     pts[] = Point3f.(x[1:n], y[1:n], z[1:n])   # updates one observable atomically
# end



# #calcuatin the sub-cubes interpolations stored in the Ogden_3D_mdbm.ncubes[i].posinterp
# interpsubcubesolution!(Ogden_3D_mdbm)
# #extracting the resutls to from the 
# path2points = extract_paths(Ogden_3D_mdbm);
# 
# ##extracting the unique points and plotting
# #puniq = unique(collect(Iterators.flatten(Iterators.flatten(path2points))));
# #scatter!(ax2, getindex.(puniq, 1), getindex.(puniq, 2), getindex.(puniq, 3), markersize=5, color=:green, label="subface - solution")
# 
# 
# #exctracing the simplexes for each ncube
# flatened_path2points = collect(Iterators.flatten(path2points))
# #eliminating the points with less than 2 points (caused by fininte precision)
# true_truflatened_path2points = flatened_path2points[length.(flatened_path2points).==3]
# 
# 
# #plotting the lines between the points
# n_faces = reshape(1:(3*length(true_truflatened_path2points)), (3, length(true_truflatened_path2points)))'
# vertices_mat = hcat(Iterators.flatten(true_truflatened_path2points)...)
# mesh!(vertices_mat, n_faces, alpha=0.5, label="subface - local simplex")




## --------------------- 3D case - 2 equtaion ------------------------


R_Ogden_λλTv0(λ::Float64, λT::Float64, v0::Float64) = R_Ogden(λ, λT, v0, α_fix)::Float64
dR_Ogden_λλTv0(λ::Float64, λT::Float64, v0::Float64) = dR_Ogden_dλT(λ, λT, v0, α_fix)::Float64
#dR_Ogden_λλTv0(λ::Float64, λT::Float64, v0::Float64)=dR_Ogden_dλT(λ, λT, v0,α_fix,dofinitediff=true)::Float64


OgdenComib(λ::Float64, λT::Float64, v0::Float64) = (R_Ogden(λ, λT, v0, α_fix), dR_Ogden_dλT(λ, λT, v0, α_fix))::Tuple{Float64,Float64}

Ogden_3D_combi_mdbm = MDBM_Problem(OgdenComib, [axis_λ, axis_λT, axis_v0])
@time solve!(Ogden_3D_combi_mdbm, 6, interpolationorder=1, checkneighbourNum=2)#number of refinements - increase it slightly to see smoother results 


#calcuatin the sub-cubes interpolations stored in the Ogden_3D_mdbm.ncubes[i].posinterp
interpsubcubesolution!(Ogden_3D_combi_mdbm)
#extracting the resutls to from the 
path2points = extract_paths(Ogden_3D_combi_mdbm);

#extracting the unique points and plotting
puniq = unique(collect(Iterators.flatten(Iterators.flatten(path2points))));
scatter!(ax2, getindex.(puniq, 1), getindex.(puniq, 2), getindex.(puniq, 3), markersize=12, color=:black, label="subface - solution")





f = Figure(size=(1400, 1000))
#show the final resolution of the grid based on the minorticks
kwargs = (; xminorticksvisible=true, xminorgridvisible=true, yminorticksvisible=true, yminorgridvisible=true)
ax1 = GLMakie.Axis(f[1, 1]; kwargs..., xlabel="λ", ylabel="v0")
scatter!(ax1, getindex.(puniq, 1), getindex.(puniq, 3), markersize=8, color=:black, label="subface - solution")
f














## --------------------- 3D case - 1 equtaion ------------------------



R_Ogden_λλTα(λ::Float64, λT::Float64, α::Float64) = R_Ogden(λ, λT, v0_fix, α)::Float64
Ogden_3D_mdbm2 = MDBM_Problem(R_Ogden_λλTα, [axis_λ, axis_λT, axis_alfa])
@time solve!(Ogden_3D_mdbm2, 6, interpolationorder=1, checkneighbourNum=0)#number of refinements - increase it slightly to see smoother results 


f = Figure(size=(1400, 1000))
ax2 = GLMakie.Axis3(f[1, 1], xlabel="λ", ylabel="λT", zlabel="α")


x, y, z = getinterpolatedsolution(Ogden_3D_mdbm2)
scatter!(ax2, x, y, z, markersize=3, label="solution")

# #calcuatin the sub-cubes interpolations stored in the Ogden_3D_mdbm2.ncubes[i].posinterp
#  interpsubcubesolution!(Ogden_3D_mdbm2)
#  #extracting the resutls to from the 
#  path2points = extract_paths(Ogden_3D_mdbm2);
#  
#  ##extracting the unique points and plotting
#  #puniq = unique(collect(Iterators.flatten(Iterators.flatten(path2points))));
#  #scatter!(ax2, getindex.(puniq, 1), getindex.(puniq, 2), getindex.(puniq, 3), markersize=5, color=:green, label="subface - solution")
#  
#  
#  
#  #exctracing the simplexes for each ncube
#  flatened_path2points = collect(Iterators.flatten(path2points))
#  #eliminating the points with less than 2 points (caused by fininte precision)
#  true_truflatened_path2points = flatened_path2points[length.(flatened_path2points).==3]
#  
#  
#  #plotting the lines between the points
#  n_faces = reshape(1:(3*length(true_truflatened_path2points)), (3, length(true_truflatened_path2points)))'
#  vertices_mat = hcat(Iterators.flatten(true_truflatened_path2points)...)
#  mesh!(vertices_mat, n_faces, alpha=0.5, label="subface - local simplex")




## --------------------- 3D case - 2 equtaion ------------------------



R_Ogden_λλTα(λ::Float64, λT::Float64, α::Float64) = R_Ogden(λ, λT, v0_fix, α)::Float64
dR_Ogden_λλTα(λ::Float64, λT::Float64, α::Float64) = dR_Ogden_dα(λ, λT, v0_fix, α)::Float64
#dR_Ogden_λλTv0(λ::Float64, λT::Float64, v0::Float64)=dR_Ogden_dλT(λ, λT, v0,α_fix,dofinitediff=true)::Float64


OgdenComib(λ::Float64, λT::Float64, α::Float64) = (R_Ogden_λλTα(λ, λT, α), dR_Ogden_λλTα(λ, λT, α))::Tuple{Float64,Float64}

Ogden_3D_combi_mdbm = MDBM_Problem(OgdenComib, [axis_λ, axis_λT, axis_alfa])
@time solve!(Ogden_3D_combi_mdbm, 6, interpolationorder=1, checkneighbourNum=2)#number of refinements - increase it slightly to see smoother results 


#calcuatin the sub-cubes interpolations stored in the Ogden_3D_mdbm.ncubes[i].posinterp
interpsubcubesolution!(Ogden_3D_combi_mdbm)
#extracting the resutls to from the 
path2points = extract_paths(Ogden_3D_combi_mdbm);

#extracting the unique points and plotting
puniq = unique(collect(Iterators.flatten(Iterators.flatten(path2points))));
scatter!(ax2, getindex.(puniq, 1), getindex.(puniq, 2), getindex.(puniq, 3), markersize=12, color=:black, label="subface - solution")





f = Figure(size=(1400, 1000))
#show the final resolution of the grid based on the minorticks
kwargs = (; xminorticksvisible=true, xminorgridvisible=true, yminorticksvisible=true, yminorgridvisible=true)
ax1 = GLMakie.Axis(f[1, 1]; kwargs..., xlabel="λ", ylabel="α")
scatter!(ax1, getindex.(puniq, 1), getindex.(puniq, 3), markersize=8, color=:black, label="subface - solution")
f
















# ===================================================================================================================================
# ===================================================================================================================================
# ===================================================================================================================================


## --------------------- 4D case - 2 equtaion ------------------------



axis_λ = MDBM.Axis(collect(LinRange(0.05, 0.5, 5)), "λ")
axis_λT = MDBM.Axis(collect(LinRange(0.05, 3, 5)), "λT")
axis_v0 = MDBM.Axis(collect(LinRange(0.25, 0.49, 5)), "v0")
axis_alfa = MDBM.Axis(collect(LinRange(-3.0, 8.0, 8)), "α")


OgdenComib4D(λ::Float64, λT::Float64, v0::Float64, α::Float64) = (R_Ogden(λ, λT, v0, α), dR_Ogden_dλT(λ, λT, v0, α))::Tuple{Float64,Float64}
#OgdenComib4D(λ::Float64, λT::Float64, v0::Float64,α::Float64)=(R_Ogden(λ, λT,v0,α ),dR_Ogden_dλT(λ, λT, v0,α,dofinitediff=true,h=1e-4))::Tuple{Float64,Float64}


Ogden_4D_combi_mdbm = MDBM_Problem(OgdenComib4D, [axis_λ, axis_λT, axis_v0, axis_alfa])
@time solve!(Ogden_4D_combi_mdbm, 6, interpolationorder=1, normp=1e3, ncubetolerance=-0.35, checkneighbourNum=0)

 checkneighbour!(Ogden_4D_combi_mdbm; interpolationorder=1, maxiteration=5, normp=20.0, ncubetolerance=0.5, doThreadprecomp=true)


f3 = Figure(size=(1400, 1000))
ax3 = GLMakie.Axis3(f3[1, 1], xlabel="λ", ylabel="v0", zlabel="α")
# # n-cube interpolation
x, y, z, w = getinterpolatedsolution(Ogden_4D_combi_mdbm)
scatter!(ax3, x, z, w, markersize=3, color=y, label="solution")
Colorbar(f3[1, 2])
# 

f3


f = Figure(size=(1400, 1000))
#show the final resolution of the grid based on the minorticks
kwargs = (; xminorticksvisible=true, xminorgridvisible=true, yminorticksvisible=true, yminorgridvisible=true)
ax1 = GLMakie.Axis(f[1, 1]; kwargs..., xlabel="λ", ylabel="v0")
scatter!(ax1, x, z, markersize=8, color=:black, label="subface - solution")


# # connecting and plotting the "mindpoints" of the n-cubes
# DT1 = connect(Ogden_4D_combi_mdbm)
# edge2plot_xyz = [reduce(hcat, [i_sol[getindex.(DT1, 1)], i_sol[getindex.(DT1, 2)], fill(NaN, length(DT1))])'[:] for i_sol in xyz_sol]
# lines!(edge2plot_xyz..., linewidth=5, label="midpoints solution connected")
# 



## --------------------- 4D case - 3 equtaion ------------------------

OgdenComib4D_extra(λ::Float64, λT::Float64, v0::Float64, α::Float64) =
    (R_Ogden(λ, λT, v0, α), dR_Ogden_dλT(λ, λT, v0, α), dR_Ogden_dα(λ, λT, v0, α))::Tuple{Float64,Float64,Float64}

# OgdenComib4D_extra(λ::Float64, λT::Float64, v0::Float64,α::Float64)=
# (R_Ogden(λ, λT,v0,α ),dR_Ogden_dα(λ, λT, v0,α,dofinitediff=true))::Tuple{Float64,Float64}

#OgdenComib4D_extra(λ::Float64, λT::Float64, v0::Float64,α::Float64)=
#(R_Ogden(λ, λT,v0,α ),dR_Ogden_dλT(λ, λT, v0,α,dofinitediff=true),dR_Ogden_dα(λ, λT, v0,α,dofinitediff=true))::Tuple{Float64,Float64,Float64}

Ogden_4D_combi_ext_mdbm = MDBM_Problem(OgdenComib4D_extra, [axis_λ, axis_λT, axis_v0, axis_alfa])
#@time solve!(Ogden_4D_combi_ext_mdbm, 3, interpolationorder=1, checkneighbourNum=2)
@time solve!(Ogden_4D_combi_ext_mdbm, 6, interpolationorder=0, normp=1e3, ncubetolerance=1.5, checkneighbourNum=0, doThreadprecomp=true)


#@time solve!(Ogden_4D_combi_ext_mdbm, 3, interpolationorder=1, checkneighbourNum=0)

# @time solve!(Ogden_4D_combi_ext_mdbm, 1, interpolationorder=1, checkneighbourNum=0) 




x_4p3, y_4p3, z_4p3, w_4p3 = getinterpolatedsolution(Ogden_4D_combi_ext_mdbm)
scatter!(ax3, x_4p3, z_4p3, w_4p3, markersize=8, color=:red, label="solution")
# 
f3

x_4p3, y_4p3, z_4p3, w_4p3 = getinterpolatedsolution(Ogden_4D_combi_ext_mdbm)
scatter!(ax1, x_4p3, z_4p3, markersize=9, color=:red, label="solution")
f



# # connecting and plotting the "mindpoints" of the n-cubes
# DT1 = connect(Ogden_4D_combi_mdbm)
# edge2plot_xyz = [reduce(hcat, [i_sol[getindex.(DT1, 1)], i_sol[getindex.(DT1, 2)], fill(NaN, length(DT1))])'[:] for i_sol in xyz_sol]
# lines!(edge2plot_xyz..., linewidth=5, label="midpoints solution connected")
# 
