
5+5
using MDBM
using SemiDiscretizationMethod


function createMillingProblem(n_spindlespeed, a_axialdeps)
    #TODO: replace these matrixes with user defined ones
    AMx = ProportionalMX(t -> @SMatrix [0.0 1.0; -δ-ε*cos(2π / T * t) -a1])
    τ1 = t -> 2π # if function is needed, the use τ1 = t->foo(t)
    BMx1 = DelayMX(τ1, t -> @SMatrix [0.0 0.0; b0 0.0])
    cVec = Additive(t -> @SVector [0.0, sin(4π / T * t)])
    LDDEProblem(AMx, [BMx1], cVec)
end;


using MDBM

method=SemiDiscretization(2,T/40);
taumax=1 #TODO: replace with the correct walue
@warn "The value of τmax is set to 1, which may not be correct. Please replace it with the correct value."
foo(n::Float64,a::Float64)::Float64 = log(spectralRadiusOfMapping(DiscreteMapping_LR(createMillingProblem(δ,ε,b0,a1,T=T),method,τmax,
    n_steps=Int((T+100eps(T))÷method.Δt)))); # No additive term calculated
using MDBM
axis=[MDBM.Axis(-1.01:0.2:5.,:δ),
    MDBM.Axis(-2.01:0.2:1.5,:b0)]

iteration=3;
@time foo(0.3,0.3)
ns=-1.01:0.1:5.0#TDOO:replace
as=-2.01:0.1:1.5#TDOO:replace

@time stab_border_points=getinterpolatedsolution(solve!(createMillingProblem(foo,[δs,b0s]),3,doThreadprecomp=false));

scatter(stab_border_points...,xlim=(-1.,5),ylim=(-2.,1.5),
    label="",    guidefontsize=14,tickfont = font(10),markersize=2,markerstrokewidth=0)