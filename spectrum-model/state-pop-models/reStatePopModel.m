function [out] = reStatePopModel(t,v,det,P,b,Om,dEz,gamL,gamD)
% t (1x2 double): time window to return state pop solutions in units of gamma^-1: t(1) and t(2) are start and end times, respectively
% v (double): mean ion velocity with units gam422/k
% det (double): laser detuning from resonance with units gam422
% p (double): initial spin-polarization of state |2>, everything else in |1>
% c (struct): useful constants in SI units
% b (struct): basis state information, see 'basis/defineBasisStates.m'
% eps (double): eps(i,j) is polarization for transition |j> to |i>
% Om (4x4 complex double): Rabi coupling
% dEz (1x4 double): dEz(k) is the Zeeman shift (with units \hbargam422) of basis state |k> 
% gamL (double): LIF-laser linewidth in units gam422
% gamD (double): decay rate into 2D_{3/2} state units gam422

R = getScatteringRate(det,v,dEz,Om,gamL,b);

% gam(i,j) is the decay rate for transition |i> to |j>
gam = zeros(b.numstates);
for k = 1:b.numstates
    gam(b.exc(k),b.gnd(k)) = b.gam(k);
end

% specify initial conditions
P0 = [(1-P)/2 (1+P)/2 0 0]';

% define rate equations - note all input frequencies in units gam
zerovec = zeros(b.numstates,1);
[dpdt] = @(t,p) defineRateEquations(t,p,R,gam,zerovec,gamD);

% solve equations using built in Matlab solver - use to check accuracy of custom rk integrator, which is faster.
% tol = 1e-6;
% opts = odeset('RelTol',tol);
% [tSol,pSol] = ode45(dpdt,[t(1) t(end)],p0,opts);
% tSol = tSol';
% pSol = pSol';

% solve equations using custom rk integrator
[tSol,pSol] = ode_rk_integrator(t,P0,dpdt);

out.det = det;
out.t = tSol;
out.p = pSol;

end