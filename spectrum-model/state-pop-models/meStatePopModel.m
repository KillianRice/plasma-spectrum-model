function [out] = meStatePopModel(t,v,det,B,P,tol,c,b,Om,dEz)
% t (1xn double): time vector for solutions in units of gamma^-1
% v (1x1 double): mean ion velocity with units \gamma/k
% det (1x1 double): laser detuning from ionic resonance
% I (1x1 double): uniform laser intensity with units W/m^2
% B (1x1 double): magnitude of magnetic field in Tesla
% theta (1x1 double): relative angle between laser polarization and ion quantization axis in radians
% P (1x1 double): electropn-spin polarization of 2S_{1/2} ground state of Sr+
% pol (string): polarization of imaging laser ('lin', 'left', 'right')
% tol (1x1 double): tolerance for differential equation solver (typical 1e-6)
% c (struct): useful constants see 'basis/defineConstants.m'
% b (struct): basis state information, see 'basis/defineBasisStates.m'
% eps (double): eps(i,j) is polarization for transition |j> to |i>
% Om0 (1x1 double): State-independent Rabi freq (units gam422)
% Om (nxn complex double): Rabi coupling
% dEz (1xn double): Zeeman shift (with units \hbar\gamma) of basis state |n> 

% This function solves the master equation for the open, four-level LIF configuration for the coherent time
% evolution of the density matrix including effects of both Zeeman and Doppler shifts and the dipole matrix
% element of a linearly-polarized laser whose polarization subtends angle theta from the local quantization
% axis. For more information, please see the OneNote entry from 03/18/20.
[H] = defineHamiltonian(v,det,B,Om,b,c);


%% Define and solve differential equations for master equation
% define master equations
inpg = b.gnd;
inpe = b.exc;
inpgam = b.gam;
inpzeromat = complex(zeros(b.numstates),zeros(b.numstates));

masterEquations = @(t,p) defineMasterEquation(t,p,H,inpg,inpe,inpgam,inpzeromat);

% construct initial conditions
p0 = complex(zeros(b.numstates),zeros(b.numstates));
p0(1,1) = (1-P)/2;
p0(2,2) = (1+P)/2;
p0 = p0(:);

p0R = real(p0);
p0I = imag(p0);
p0In = [p0R;p0I];

% define options for ode solver
opts = odeset('RelTol',tol);

% solve master equations using ode45
[tSol,pSol] = ode45(masterEquations,t,p0In,opts);


%% Output Parameters

out.det = det;
out.t = tSol;
out.p = pSol;

out.det = det;
out.t = tSol;
out.p = pSol(:,[1,6,11,16])';


end