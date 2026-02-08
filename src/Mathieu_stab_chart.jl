
5+5
using MDBM
using SemiDiscretizationMethod


function createMathieuProblem(δ, ε, b0, a1; T=2π)
    AMx = ProportionalMX(t -> @SMatrix [0.0 1.0; -δ-ε*cos(2π / T * t) -a1])
    τ1 = t -> 2π # if function is needed, the use τ1 = t->foo(t)
    BMx1 = DelayMX(τ1, t -> @SMatrix [0.0 0.0; b0 0.0])
    cVec = Additive(t -> @SVector [0.0, sin(4π / T * t)])
    LDDEProblem(AMx, [BMx1], cVec)
end;

τmax=2π # the largest τ of the system
T=2π #Principle period of the system (sin(t)=sin(t+T)) 
mathieu_lddep=createMathieuProblem(3.,2.,-0.15,0.1,T=T); # LDDE problem for Hayes equation
method=SemiDiscretization(1,0.01) # 1st order semi discretization with Δt=0.01
# if P = τmax, then n_steps is automatically calculated
mapping=DiscreteMapping_LR(mathieu_lddep,method,τmax,
    n_steps=Int((T+100eps(T))÷method.Δt),calculate_additive=true); #The discrete mapping of the system

@time spectralRadiusOfMapping(mapping); # spectral radius ρ of the mapping matrix (ρ>1 unstable, ρ<1 stable)
# spectralRadiusOfMapping(mapping) = 0.5131596340374617

@time fp=fixPointOfMapping(mapping); # stationary solution of the hayes equation (equilibrium position)

# fixPointOfMapping plotted
using Plots
gr();
using LaTeXStrings
n_steps=Int((T+100eps(T))÷method.Δt)
t_plot=0.0:method.Δt:((size(fp,1)/2-1)*method.Δt)
plot(t_plot,fp[1:2:end],
    xlabel=L"-s",title=L"t \in [nT,(n+1)T],\quad n \to \infty",guidefontsize=14,linewidth=3,
    label=L"x(t-s)",legendfontsize=11,tickfont = font(10))
plot!(t_plot,fp[2:2:end],
    xlabel=L"-s",linewidth=3,
    label=L"\dot{x}(t-s)")
plot!(t_plot,sin.(2 .* t_plot),linewidth=3,label=L"\sin(2t)")



using MDBM
a1=0.1;
ε=1;
τmax=2π;
T=1π;
method=SemiDiscretization(2,T/40);

foo(δ::Float64,b0::Float64)::Float64 = log(spectralRadiusOfMapping(DiscreteMapping_LR(createMathieuProblem(δ,ε,b0,a1,T=T),method,τmax,
    n_steps=Int((T+100eps(T))÷method.Δt)))); # No additive term calculated
using MDBM
axis=[MDBM.Axis(-1.01:0.2:5.,:δ),
    MDBM.Axis(-2.01:0.2:1.5,:b0)]

iteration=3;
@time foo(0.3,0.3)
δs=-1.01:0.1:5.0
b0s=-2.01:0.1:1.5
@show length(δs),length(b0s)
@show length(δs)*length(b0s)
@time mus=[foo(δ,b0) for δ in δs, b0 in b0s];
contour(δs,b0s,mus',levels=[0.0],xlabel=L"\delta",ylabel=L"b_0",
    title="Stability border of the delay Mathieu equation",guidefontsize=14,tickfont = font(10),linewidth=3)

    
δs=-1.01:0.5:5.0
b0s=-2.01:0.5:1.5
iteration=3;
@time stab_border_points=getinterpolatedsolution(solve!(MDBM_Problem(foo,[δs,b0s]),3,doThreadprecomp=false));

scatter(stab_border_points...,xlim=(-1.,5),ylim=(-2.,1.5),
    label="",title="Stability border of the delay Mathieu equation",xlabel=L"\delta",ylabel=L"b_0",
    guidefontsize=14,tickfont = font(10),markersize=2,markerstrokewidth=0)