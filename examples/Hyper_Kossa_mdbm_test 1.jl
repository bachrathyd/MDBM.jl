# Attila Kossa  -- Hyper stability example - using MDBM to find the solutions

5 + 5
using Revise
using MDBM
using GLMakie

using LinearAlgebra

using ForwardDiff
GLMakie.closeall()
GLMakie.activate!(; title="2 parameters, codimension 1")
# -----------------------------

# Unicode version
function R_Hyper(X, Y, B, P)
    return (2.0 / 3.0) .*
           ((X .* Y .^ 2.0) .^ (-5.0 / 3.0) .* ((-1.0) .* X .^ 2.0 + Y .^ 2.0) + (-2.0) .* (1.0 + B) .* (1.0 + P) .* ((-1.0) + 2.0 .* P) .^ (-1.0) .*
           ((-1.0) + X .* Y .^ 2.0) + B .* (X .* Y .^ 2.0) .^ (-7.0 / 3.0) .* ((-1.0) .* X .^ 2.0 .* Y .^ 2.0 + Y .^ 4.0))
end

"""
Wrapper for R_Hyper that uses a parameter β such that B = tan(β).
This allows B to span from -Inf to +Inf.
"""
function R_Hyper_beta(X, Y, β, P)
    B = tan(β)
    return R_Hyper(X, Y, B, P)
end





"""
Calculates the derivative of R_Hyper with respect to Y.

By default, it uses automatic differentiation (ForwardDiff).
If `dofinitediff=true` is passed, it uses a central finite difference
method with step size `h`.
"""
function dR_Hyper_dY(X::Float64,
    Y::Float64,
    B::Float64,
    P::Float64;  # Semicolon separates positional from keyword args
    dofinitediff::Bool=false,
    h::Float64=1e-2)::Float64

    if dofinitediff
        # --- Central Finite Difference Method ---
        # f'(x) ≈ (f(x+h) - f(x-h)) / (2h)        
        f_plus = R_Hyper(X, Y + h, B, P)
        f_minus = R_Hyper(X, Y - h, B, P)

        return (f_plus - f_minus) / (2.0 * h)

    else
        # --- Automatic Differentiation (ForwardDiff) Method ---
        # Original code
        # We create an anonymous function that only takes the parameter of interest
        return ForwardDiff.derivative(lamdaT -> R_Hyper(X, lamdaT, B, P), Y)

    end
end




"""
Calculates the derivative of R_Hyper with respect to P.

By default, it uses automatic differentiation (ForwardDiff).
If `dofinitediff=true` is passed, it uses a central finite difference
method with step size `h`.
"""
function dR_Hyper_dP(X::Float64,
    Y::Float64,
    B::Float64,
    P::Float64;  # Semicolon separates positional from keyword args
    dofinitediff::Bool=false,
    h::Float64=1e-0)::Float64

    if dofinitediff
        # --- Central Finite Difference Method ---
        # f'(x) ≈ (f(x+h) - f(x-h)) / (2h)

        f_plus = R_Hyper(X, Y, B, P + h)
        f_minus = R_Hyper(X, Y, B, P - h)

        return (f_plus - f_minus) / (2.0 * h)

    else
        # --- Automatic Differentiation (ForwardDiff) Method ---
        # Original code
        # We create an anonymous function that only takes the parameter of interest
        f = P_par -> R_Hyper(X, Y, B, P_par)

        return ForwardDiff.derivative(f, P)
    end
end

const P_fix = 0.41 # Fixed value for P
const β_fix = atan(-0.09) # Fixed value for β, corresponding to B = -0.09


R_Hyper_XY(X::Float64, Y::Float64) = R_Hyper_beta(X, Y, β_fix, P_fix)::Float64
dR_Hyper_dlambda_XY(X::Float64, Y::Float64) = dR_Hyper_dY(X, Y, tan(β_fix), P_fix)::Float64
#dR_Hyper_dlambda_XY(X::Float64, Y::Float64)=dR_Hyper_dY(X, Y, B_fix,P_fix,dofinitediff=true)::Float64


axis_X = MDBM.Axis((0.01:0.101:1.2) .^ 2, "X")
axis_Y = MDBM.Axis((0.01:0.201:2.0) .^ 2, "Y")

Hyper_mdbm = MDBM_Problem(R_Hyper_XY, [axis_X, axis_Y])
@time solve!(Hyper_mdbm, 4)#number of refinements - increase it slightly to see smoother results 
Hyper_deiv_mdbm = MDBM_Problem(dR_Hyper_dlambda_XY, [axis_X, axis_Y])
@time solve!(Hyper_deiv_mdbm, 4)#number of refinements - increase it slightly to see smoother results 



f = Figure()
#show the final resolution of the grid based on the minorticks
kwargs = (; xminorticksvisible=true, xminorgridvisible=true, yminorticksvisible=true, yminorgridvisible=true)
ax1 = GLMakie.Axis(f[1, 1]; xminorticks=Hyper_mdbm.axes[1].ticks, yminorticks=Hyper_mdbm.axes[2].ticks, kwargs..., xlabel="X", ylabel="Y")

# # n-cube interpolation
xy_sol = getinterpolatedsolution(Hyper_mdbm)
# scatter!(xy_sol..., markersize = 15, color = :red,marker ='x',strokewidth=3,label = "solution")
# 
# # show the points where the function is evaluated
# xy_val = getevaluatedpoints(Hyper_mdbm)
# fval=getevaluatedfunctionvalues(Hyper_mdbm)
# scatter!(xy_val...,color=sign.(fval),label = "evaluated")

# connecting and plotting the "mindpoints" of the n-cubes
DT1 = connect(Hyper_mdbm)
edge2plot_xyz = [reduce(hcat, [i_sol[getindex.(DT1, 1)], i_sol[getindex.(DT1, 2)], fill(NaN, length(DT1))])'[:] for i_sol in xy_sol]
lines!(edge2plot_xyz..., linewidth=5, label="midpoints solution connected")

# connecting and plotting the "mindpoints" of the n-cubes
xy_sol = getinterpolatedsolution(Hyper_deiv_mdbm)
DT1 = connect(Hyper_deiv_mdbm)
edge2plot_xyz = [reduce(hcat, [i_sol[getindex.(DT1, 1)], i_sol[getindex.(DT1, 2)], fill(NaN, length(DT1))])'[:] for i_sol in xy_sol]
lines!(edge2plot_xyz..., linewidth=5, label="midpoints solution connected")

f
## --------------------- 3D case - 1 equtaion ------------------------



const P_fix = 0.3


axis_X = MDBM.Axis(LinRange(0.25, 1.2, 10) .^ 2, "X")
axis_Y = MDBM.Axis(LinRange(0.25, 2.0, 10) .^ 2, "Y")
axis_β = MDBM.Axis(LinRange(atan(-1.401), atan(1.02), 6), "β")
axis_β = MDBM.Axis(LinRange(atan(-100.401), atan(100.02), 6), "β")
axis_β = MDBM.Axis(LinRange(-2.0, 4.0, 16), "β")


R_Hyper_XYβ(X::Float64, Y::Float64, β::Float64) = R_Hyper_beta(X, Y, β, P_fix)::Float64

Hyper_3D_mdbm = MDBM_Problem(R_Hyper_XYβ, [axis_X, axis_Y, axis_β])
@time solve!(Hyper_3D_mdbm, 3, interpolationorder=1,checkneighbourNum=0)#number of refinements - increase it slightly to see smoother results 


f = Figure()
ax2 = GLMakie.Axis3(f[1, 1], xlabel="X", ylabel="Y", zlabel="β")


x,y,z = getinterpolatedsolution(Hyper_3D_mdbm)
scatter!(ax2,x,y,z, markersize = 3,label = "solution")
# #calcuatin the sub-cubes interpolations stored in the Hyper_3D_mdbm.ncubes[i].posinterp
# interpsubcubesolution!(Hyper_3D_mdbm)
# #extracting the resutls to from the 
# path2points = extract_paths(Hyper_3D_mdbm);
# 
# ##extracting the unique points and plotting
# #puniq = unique(collect(Iterators.flatten(Iterators.flatten(path2points))));
# #scatter!(ax2, getindex.(puniq, 1), getindex.(puniq, 2), getindex.(puniq, 3), markersize=5, color=:green, label="subface - solution")
# 
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


axis_X = MDBM.Axis(LinRange(0.25, 1.2, 10) .^ 2, "X")
axis_Y = MDBM.Axis(LinRange(0.25, 2.0, 10) .^ 2, "Y")
axis_β = MDBM.Axis(LinRange(atan(-1.401), atan(1.02), 6), "β")



R_Hyper_XYβ(X::Float64, Y::Float64, β::Float64) = R_Hyper_beta(X, Y, β, P_fix)::Float64
dR_Hyper_XYβ(X::Float64, Y::Float64, β::Float64) = dR_Hyper_dY(X, Y, tan(β), P_fix)::Float64
#dR_Hyper_XYB(X::Float64, Y::Float64, B::Float64)=dR_Hyper_dY(X, Y, B,P_fix,dofinitediff=true)::Float64


HyperComib(X::Float64, Y::Float64, β::Float64) = (R_Hyper_beta(X, Y, β, P_fix), dR_Hyper_dY(X, Y, tan(β), P_fix))::Tuple{Float64,Float64}

Hyper_3D_combi_mdbm = MDBM_Problem(HyperComib, [axis_X, axis_Y, axis_β])
@time solve!(Hyper_3D_combi_mdbm, 4, interpolationorder=1)#number of refinements - increase it slightly to see smoother results 


#calcuatin the sub-cubes interpolations stored in the Hyper_3D_mdbm.ncubes[i].posinterp
interpsubcubesolution!(Hyper_3D_combi_mdbm)
#extracting the resutls to from the 
path2points = extract_paths(Hyper_3D_combi_mdbm);

#extracting the unique points and plotting
puniq = unique(collect(Iterators.flatten(Iterators.flatten(path2points))));
scatter!(ax2, getindex.(puniq, 1), getindex.(puniq, 2), getindex.(puniq, 3), markersize=12, color=:black, label="subface - solution")





f5 = Figure()
#show the final resolution of the grid based on the minorticks
kwargs = (; xminorticksvisible=true, xminorgridvisible=true, yminorticksvisible=true, yminorgridvisible=true)
ax1 = GLMakie.Axis(f5[1, 1]; kwargs..., xlabel="X", ylabel="β")
scatter!(ax1, getindex.(puniq, 1), getindex.(puniq, 3), markersize=8, color=:black, label="subface - solution")
f5




# # ===================================================================================================================================
# # ===================================================================================================================================
# # ===================================================================================================================================
# 
# 
# ## --------------------- 4D case - 2 equtaion ------------------------
# 
# 
# 
# axis_X = MDBM.Axis(LinRange(0.25, 1.2, 10) .^ 2, "X")
# axis_Y = MDBM.Axis(LinRange(0.25, 2.0, 10) .^ 2, "Y")
# axis_B = MDBM.Axis(LinRange(-1.401, 1.02, 6), "B")
# axis_alfa = MDBM.Axis(0.401:0.01:0.49, "P")
# 
# HyperComib4D(X::Float64, Y::Float64, B::Float64, P::Float64) = (R_Hyper(X, Y, B, P), dR_Hyper_dY(X, Y, B, P))::Tuple{Float64,Float64}
# #HyperComib4D(X::Float64, Y::Float64, B::Float64,P::Float64)=(R_Hyper(X, Y,B,P ),dR_Hyper_dY(X, Y, B,P,dofinitediff=true,h=1e-4))::Tuple{Float64,Float64}
# 
# 
# Hyper_4D_combi_mdbm = MDBM_Problem(HyperComib4D, [axis_X, axis_Y, axis_B, axis_alfa])
# @time solve!(Hyper_4D_combi_mdbm, 3, interpolationorder=1, checkneighbourNum=0)#number of refinements - increase it slightly to see smoother results 
# 
# 
# 
# f3 = Figure()
# ax3 = GLMakie.Axis3(f3[1, 1], xlabel="X", ylabel="β", zlabel="P")
# # # n-cube interpolation
# x, y, z, w = getinterpolatedsolution(Hyper_4D_combi_mdbm)
# scatter!(ax3, x, z, w, markersize=5, color=y, label="solution")
# Colorbar(f3[1, 2])
# # 
# 
# 
# 
# f = Figure()
# #show the final resolution of the grid based on the minorticks
# kwargs = (; xminorticksvisible=true, xminorgridvisible=true, yminorticksvisible=true, yminorgridvisible=true)
# ax1 = GLMakie.Axis(f[1, 1]; kwargs..., xlabel="X", ylabel="β")
# scatter!(ax1, x, z, markersize=8, color=:black, label="subface - solution")
# 
# 
# # # connecting and plotting the "mindpoints" of the n-cubes
# # DT1 = connect(Hyper_4D_combi_mdbm)
# # edge2plot_xyz = [reduce(hcat, [i_sol[getindex.(DT1, 1)], i_sol[getindex.(DT1, 2)], fill(NaN, length(DT1))])'[:] for i_sol in xyz_sol]
# # lines!(edge2plot_xyz..., linewidth=5, label="midpoints solution connected")
# # 
# 
# 
# 
# ## --------------------- 4D case - 3 equtaion ------------------------
# 
# 
# axis_X = MDBM.Axis(LinRange(0.25, 1.2, 10) .^ 2, "X")
# axis_Y = MDBM.Axis(LinRange(0.25, 2.0, 10) .^ 2, "Y")
# axis_B = MDBM.Axis(LinRange(-1.401, 1.02, 6), "B")
# axis_alfa = MDBM.Axis(0.401:0.01:0.49, "P")
# 
# HyperComib4D_extra(X::Float64, Y::Float64, B::Float64, P::Float64) =
#     (R_Hyper(X, Y, B, P), dR_Hyper_dY(X, Y, B, P), dR_Hyper_dP(X, Y, B, P, dofinitediff=false))::Tuple{Float64,Float64,Float64}
# 
# # HyperComib4D_extra(X::Float64, Y::Float64, B::Float64,P::Float64)=
# # (R_Hyper(X, Y,B,P ),dR_Hyper_dP(X, Y, B,P,dofinitediff=true))::Tuple{Float64,Float64}
# 
# #HyperComib4D_extra(X::Float64, Y::Float64, B::Float64,P::Float64)=
# #(R_Hyper(X, Y,B,P ),dR_Hyper_dY(X, Y, B,P,dofinitediff=true),dR_Hyper_dP(X, Y, B,P,dofinitediff=true))::Tuple{Float64,Float64,Float64}
# 
# Hyper_4D_combi_ext_mdbm = MDBM_Problem(HyperComib4D_extra, [axis_X, axis_Y, axis_B, axis_alfa])
# #@time solve!(Hyper_4D_combi_ext_mdbm, 3, interpolationorder=1, checkneighbourNum=2)
# @time solve!(Hyper_4D_combi_ext_mdbm, 3, interpolationorder=1, normp=20, ncubetolerance=1.5, checkneighbourNum=0, doThreadprecomp=true)
# 
# 
# 
# # HyperComib4D_extra(X::Float64, Y::Float64, B::Float64,P::Float64)=
# # (R_Hyper(X, Y,B,P ),dR_Hyper_dP(X, Y, B,P))::Tuple{Float64,Float64}
# # Hyper_4D_combi_ext_mdbm = MDBM_Problem(HyperComib4D_extra, [axis_X,axis_Y,axis_B,axis_alfa])
# # @time solve!(Hyper_4D_combi_ext_mdbm,4, interpolationorder=1, normp=20, ncubetolerance=0.8, checkneighbourNum=0, doThreadprecomp=true)
# 
# 
# 
# x_4p3, y_4p3, z_4p3, w_4p3 = getinterpolatedsolution(Hyper_4D_combi_ext_mdbm)
# scatter!(ax3, x_4p3, z_4p3, w_4p3, markersize=10, color=:black, label="solution")
# # 
# f3
# 
# x_4p3, y_4p3, z_4p3, w_4p3 = getinterpolatedsolution(Hyper_4D_combi_ext_mdbm)
# scatter!(ax1, x_4p3, z_4p3, markersize=9, color=:red, label="solution")
# f
# 
# f3
# 
# # # connecting and plotting the "mindpoints" of the n-cubes
# # DT1 = connect(Hyper_4D_combi_mdbm)
# # edge2plot_xyz = [reduce(hcat, [i_sol[getindex.(DT1, 1)], i_sol[getindex.(DT1, 2)], fill(NaN, length(DT1))])'[:] for i_sol in xyz_sol]
# # lines!(edge2plot_xyz..., linewidth=5, label="midpoints solution connected")
# # 
# 