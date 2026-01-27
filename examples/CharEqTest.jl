## MDBM for balancing an e-scooeter
# Linearized equations of motion, control with steering with multiple delays

# import Pkg; Pkg.add("GLMakie")
# import Pkg; Pkg.add("LaTeXStrings")

using MDBM                        # Multi‐Dimensional Bisection Method 
using GLMakie                     # High‐performance plotting             
using LinearAlgebra               # for complex operations (norm, etc.)
using LaTeXStrings

mch=9.5
# fixed parameters for e-scooter (Xiaomi M365)
mch = 9.5;
mch = 9.5;
mfork = 2.797;
mwheel1 = 2.894;
mwheel2 = 1.136;
p = 0.829;
b = 0.2418;
e = 0.0281;
epsilon = deg2rad(14.2);
R = 0.111;
xCGch = 0.4314;
xCGfork = 0.01; 
zCGch = 0.0782; 
zCGfork = 0.271586;
Aw1 = 0.01134;
Bw1 = 0.00628;
Aw2 = 0.00536;
Bw2 = 0.0029;
Ach = 0.1014;
Bch = 0.4374;
Cch = 0.3561;
Dch = -0.1074;
Afork = 0.341171;
Bfork = 0.4018;
Cfork = 0.0840292;
Dfork = 0.151596;
v = 0.0;
g = 9.81;
Kpdelta = 10.0; # optimal proportional gain for the parameter setup (based on MSD article)
Kddelta = -5.0; # optimal derivative gain for the parameter setup (based on MSD article)
tau1 = 1e-10; # for steering (delta)
tau2 = 0.01;  # for lean (phi)

const parlist=(mch,mfork,mwheel1,mwheel2,p,b,e,epsilon,R,xCGch,xCGfork,zCGch,zCGfork,Aw1,Bw1,Aw2,Bw2,Ach,Bch,Cch,Dch,Afork,Bfork,Cfork,Dfork,v,g,Kpdelta,Kddelta,tau1,tau2)
##
function char_fun_scooter_hierarchical_steering_multipledelays(Kpphi, Kdphi, ω; tau2=0.01, par=par)
    
mch,mfork,mwheel1,mwheel2,p,b,e,epsilon,R,xCGch,xCGfork,zCGch,zCGfork,Aw1,Bw1,Aw2,Bw2,Ach,Bch,Cch,Dch,Afork,Bfork,Cfork,Dfork,v,g,Kpdelta,Kddelta,tau1,tau2_fake =par

	# mass matrix
	M11 = (8*Ach+4*Afork+8*Aw1+8*Aw2+4*Cfork+4*b^2*mfork+e^2*mfork+8*mch*R^2+3*mfork*R^2+8*mwheel1*R^2+8*mwheel2*R^2+4*mfork*xCGfork^2+16*mch*R*zCGch+8*mch*zCGch^2+8*b*mfork*zCGfork+4*mfork*zCGfork^2+4*mfork*(3*b*R+e*xCGfork+3*R*zCGfork)*cos(epsilon)+4*(Afork-Cfork+mfork*(b^2+R^2-xCGfork^2+2*b*zCGfork+zCGfork^2))*cos(2*epsilon)+4*b*mfork*R*cos(3*epsilon)-4*e*mfork*xCGfork*cos(3*epsilon)+4*mfork*R*zCGfork*cos(3*epsilon)-e^2*mfork*cos(4*epsilon)+mfork*R^2*cos(4*epsilon)+4*b*e*mfork*sin(epsilon)+4*mfork*R*xCGfork*sin(epsilon)+4*e*mfork*zCGfork*sin(epsilon)-8*Dfork*sin(2*epsilon)+4*e*mfork*R*sin(2*epsilon)+8*b*mfork*xCGfork*sin(2*epsilon)+8*mfork*xCGfork*zCGfork*sin(2*epsilon)+4*b*e*mfork*sin(3*epsilon)+4*mfork*R*xCGfork*sin(3*epsilon)+4*e*mfork*zCGfork*sin(3*epsilon)+2*e*mfork*R*sin(4*epsilon))/8;
	M12 = (cos(epsilon)^3*(2*Dfork*(e+p)+2*Dch*e*1/cos(epsilon)^2-2*e*mch*xCGch*(R+zCGch)*1/cos(epsilon)^2-2*Afork*p*tan(epsilon)-2*Aw1*p*tan(epsilon)+2*Afork*(e+p)*tan(epsilon)-2*Cfork*(e+p)*tan(epsilon)-2*mwheel1*p*R^2*1/cos(epsilon)^2*tan(epsilon)+4*Dfork*p*tan(epsilon)^2-2*Dfork*(e+p)*tan(epsilon)^2-2*Aw1*p*tan(epsilon)^3-2*Cfork*p*tan(epsilon)^3-mfork*(2*(e+p)*xCGfork+e*(e+2*p+e*cos(2*epsilon))*1/cos(epsilon)-2*e*R*sin(epsilon)-2*e*(b+zCGfork)*tan(epsilon)+2*p*xCGfork*tan(epsilon)^2)*(b+zCGfork+R*1/cos(epsilon)+xCGfork*tan(epsilon)+sin(epsilon)*(e-R*tan(epsilon)))))/(2*p);
	M21 = M12;
	M22 = (4*Aw2*e^2*cos(epsilon)^2+4*Cch*e^2*cos(epsilon)^2+4*e^2*mch*xCGch^2*cos(epsilon)^2+4*Aw1*(e+p)*cos(epsilon)^2*(p+e*cos(epsilon)^2)+4*Cfork*(e+p)*cos(epsilon)^2*(p+e*cos(epsilon)^2)+4*Dfork*e*(e+p)*cos(epsilon)^3*sin(epsilon)+4*mwheel1*p^2*R^2*sin(epsilon)^2+2*Aw1*p*(e+2*p+e*cos(2*epsilon))*sin(epsilon)^2+2*Cfork*p*(e+2*p+e*cos(2*epsilon))*sin(epsilon)^2+4*Dfork*e*p*cos(epsilon)*sin(epsilon)^3-Dfork*p*(e+2*p+e*cos(2*epsilon))*sin(2*epsilon)+Dfork*(e+p)*(e+2*p+e*cos(2*epsilon))*sin(2*epsilon)-Afork*e*p*sin(2*epsilon)^2-Aw1*e*p*sin(2*epsilon)^2+Afork*e*(e+p)*sin(2*epsilon)^2+Aw1*e*(e+p)*sin(2*epsilon)^2+mfork*(2*(e+p)*xCGfork*cos(epsilon)^2+2*p*xCGfork*sin(epsilon)^2+e*cos(epsilon)*(e+2*p+e*cos(2*epsilon)-2*(b+zCGfork)*sin(epsilon)-R*sin(2*epsilon)))^2)/(4*p^2);
	MM = [M11 M12; M21 M22]

	# stiffness matrix
	K11 = -(g*(2*mch*R+mfork*R+2*mwheel1*R+2*mwheel2*R+2*mch*zCGch+2*mfork*(b+zCGfork)*cos(epsilon)+mfork*R*cos(2*epsilon)+2*mfork*xCGfork*sin(epsilon)+e*mfork*sin(2*epsilon)))/2;
	K12 = -((-3*e^2*g*mfork*R-4*e*g*R*(mfork*p+mch*xCGch))*cos(epsilon)+R*(-2*e*g*mfork*xCGfork-4*g*mfork*p*xCGfork-2*e*g*mfork*xCGfork*cos(2*epsilon)-e^2*g*mfork*cos(3*epsilon)+e*g*mfork*R*sin(epsilon)-4*g*mwheel1*p*R*sin(epsilon)+2*b*e*g*mfork*sin(2*epsilon)+2*e*g*mfork*zCGfork*sin(2*epsilon)+e*g*mfork*R*sin(3*epsilon)))/(4*p*R);
	K21 = K12;
	K22 = -(-(e*g*mfork*p*R^2)+4*g*mwheel1*p^2*R^2-2*mfork*R*(b*e*g*p+e*g*p*zCGfork)*cos(epsilon)-4*g*mwheel1*p^2*R^2*cos(2*epsilon)+2*b*e*g*mfork*p*R*cos(3*epsilon)+2*e*g*mfork*p*R*zCGfork*cos(3*epsilon)+e*g*mfork*p*R^2*cos(4*epsilon)+2*e*g*mfork*p*R*xCGfork*sin(epsilon)+8*g*mfork*p^2*R*xCGfork*sin(epsilon)+2*e^2*g*mfork*p*R*sin(2*epsilon)+4*e*g*mfork*p^2*R*sin(2*epsilon)+4*e*g*mch*p*R*xCGch*sin(2*epsilon)+2*e*g*mfork*p*R*xCGfork*sin(3*epsilon)+e^2*g*mfork*p*R*sin(4*epsilon))/(8*p^2*R);
	K0 = [K11 K12; K21 K22];

    # control matrices
    Pphi = [0 0; -Kpphi*Kpdelta 0];
    Pdelta = [0 0; 0 -Kpdelta];
    Dphi = [0 0; -Kdphi*Kpdelta 0];
    Ddelta = [0 0; 0 -Kddelta];
    
    lambda = im * ω;
    charmat = MM.*lambda^2 + K0 -(Pdelta + Ddelta.*lambda).*exp(-lambda*tau1)-(Pphi + Dphi.*lambda).*exp(-lambda*tau2);
    chareq = det(charmat);

    char_eq_real = real(chareq);
    char_eq_imag = imag(chareq);
    return char_eq_real, char_eq_imag

end
#function test
char_fun_scooter_hierarchical_steering_multipledelays(-500.0,-20.0, 15.0,  par=parlist)
@code_warntype  char_fun_scooter_hierarchical_steering_multipledelays(-500.0,-20.0, 15.0,  par=parlist)

##
# Define a small wrapper for MDBM that only uses (Kpphi,Kdphi,ω):
function char_fun3d(Kpphi, Kdphi, ω )
    #Kpphi, Kdphi, ω , tau2 = x
    return char_fun_scooter_hierarchical_steering_multipledelays(Kpphi, Kdphi, ω,  par=parlist)
end
@code_warntype char_fun3d(-500.0,-20.0, 15.0)
# Choose parameter ranges
Kpphi_min, Kpphi_max = -1100.0, 100.0
Kdphi_min, Kdphi_max = -80.0, 0.0
ω_min, ω_max = 0.0001, 150.0       

# Build coarse MDBM axes in (Kpphi, Kdphi, ω):
axis_Kpphi=LinRange(Kpphi_min, Kpphi_max, 11)
axis_Kdphi = LinRange(Kdphi_min, Kdphi_max, 11)
axis_ω = LinRange(ω_min, ω_max, 11)

# linRange (olyan mint MATLAB linspace)


axes3d = [axis_Kpphi, axis_Kdphi, axis_ω]
mdbm3d = MDBM_Problem(char_fun3d, axes3d)

# Perform a few refinements (e.g. 4 or 5) to resolve curves:
refine_iters3d = 5

#@profview  solve!(mdbm3d, 2,verbosity=1)
@time solve!(mdbm3d, refine_iters3d,verbosity=1)
println("→ Finished MDBM.")


# PLOTTING CODIM‐2 CURVES in (Kpphi,Kdphi), colored by ω
GLMakie.activate!()
import GLMakie: scatter, Axis, Figure

figure1 = Figure(size=(900, 600))
xyzr_sol = getinterpolatedsolution(mdbm3d)
ax = Axis(figure1[1, 1])
scatter!(ax, xyzr_sol[1], xyzr_sol[2], color=xyzr_sol[3], markersize=6)
ax.xlabel = L"K_{\mathrm{p}\varphi}^{\mathrm s}\, [1]"
ax.ylabel = L"K_{\mathrm{d}\varphi}^{\mathrm s}\, [\mathrm{s}]"
ax.title = L"Stabilitási\ térkép"
ax.xticklabelsize = 12
ax.yticklabelsize = 12

display(figure1)

save("stabplot.png", figure1)

## 4D -------------------------------------------------------
# Define a small wrapper for MDBM that only uses (Kpphi,Kdphi,ω):
function char_fun4d(Kpphi, Kdphi, ω ,tau2_in)
    #Kpphi, Kdphi, ω , tau2 = x
    return char_fun_scooter_hierarchical_steering_multipledelays(Kpphi, Kdphi, ω, tau2=tau2_in, par=parlist)
end
char_fun4d(-500.0,-20.0, 15.0, 0.1)
tau2_min, tau2_max = 0.05, 0.2

axis_Kpphi=LinRange(Kpphi_min, Kpphi_max, 6)
axis_Kdphi = LinRange(Kdphi_min, Kdphi_max, 6)
axis_ω = LinRange(ω_min, ω_max, 6)

axis_tau2 = LinRange(tau2_min, tau2_max, 6)

mdbm4d = MDBM_Problem(char_fun4d, [axis_Kpphi, axis_Kdphi, axis_ω,axis_tau2])

@time interpolate!(mdbm4d)
@profview refine!(mdbm4d)

solve!(mdbm4d, 2,verbosity=1,checkneighbourNum=0)
figure2 = Figure(size=(1000, 700))
ax = Axis3(figure2[1, 1],
    xlabel="Kpphi",
    ylabel="Kdphi",
    zlabel="tau2",  
    title="",
    aspect=:equal,    # use data‐ratio so the box matches the ranges
    viewmode=:fit      # fit the resulting cuboid into the frame
)
#ax1 = Axis(figure1[1, 1], xlabel="V (dimensionless)", ylabel="L (dimensionless)", title="Figure 4(a) from Takács et al. (2009), Σ=1.8, ζ=0")

xyzr_sol = getinterpolatedsolution(mdbm4d)
scatter!(ax, xyzr_sol[[1, 2, 3]]..., markersize=6, color=xyzr_sol[4])