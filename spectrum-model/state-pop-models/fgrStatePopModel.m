function [out] = fgrStatePopModel(t,v,det,P,b,Om,dEz,gamL,gamD)
% t (1xn double): time vector for solutions in units of gamma^-1
% v (1x1 double): mean ion velocity with units \gamma/k
% det (1x1 double): laser detuning from ionic resonance
% p0 (1x1 double): initial spin-polarization of state |2>, everything else in |1>
% b (struct): basis state information, see 'basis/defineBasisStates.m'
% eps (double): eps(i,j) is polarization for transition |j> to |i>
% Om (nxn complex double): Rabi coupling
% dEz (1xn double): Zeeman shift (with units \hbar\gamma) of basis state |n> 
% gamL (double): LIF-laser linewidth in units gam422
% gamD (double): decay rate into 2D_{3/2} state units gam422

% This function returns the equilibrium excited state populations using the Fermi's Golden Rule
% (FGR) approximation, which assumes the Rabi frequency is sufficiently small to not significantly
% alter the ion state populations. Under this assumption we may assume optical pumping among ground
% states does not occur. T

[R] = getScatteringRate(det,v,dEz,Om,gamL,b);

% p(i,k) is the fraction of ions that occupy state |k> at time t(i)
% the ground state populations (|1> and |2>) are specified by the spin polarization, 'p0' and do not change with time.
% the excited state populations (|3> and |4>) are calculated using a power-broadened scattering rate.
tvec = linspace(t(1),t(2),round((t(2)-t(1))/25)+2);
p = zeros(length(tvec),b.numstates);
p(1,1) = (1-P)/2;
p(1,2) = (1+P)/2;
Rsat = R./(1+gamD+2.*R);
Rout = sum(Rsat,1);
Rtot = sum(Rsat,'all');
for k = 1:length(b.exc)
    exc = b.exc(k);
    gnd = b.gnd(k);
    p(1,exc) = p(1,exc) + Rsat(exc,gnd).*p(1,gnd);
end

for i = 2:length(tvec)
    p_e_tot = sum(p(i-1,b.upInd));
    p_g_tot = sum(p(i-1,b.lowInd));
    for gnd = b.lowInd
        p(i,gnd) = p(i-1,gnd) - Rout(gnd)/Rtot*p_e_tot*gamD*(tvec(i)-tvec(i-1))*2*p(i-1,gnd)/p_g_tot;
    end
    
    for k = 1:length(b.exc)
        exc = b.exc(k);
        gnd = b.gnd(k);
        p(i,exc) = p(i,exc) + Rsat(exc,gnd).*p(i,gnd);
    end
end

out.det = det;
out.t = tvec;
out.p = p';

end