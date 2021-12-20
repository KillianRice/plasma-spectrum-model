function [H] = defineHamiltonian(v,det,B,Om,b,c)
%% Input Parameters
% v (1x1 double): mean ion velocity with units \gamma/k
% det (1x1 double): laser detuning from resonance with units \gamma
% B (1x1 double): magnitude of local magnetic field with units Tesla
% Om (nxn complex double): Rabi frequency of LIF laser, with j-dep. factors removed
% tensor notation - n is the number of basis states
% b (struct): basis state information
% c (struct): important constants

zeromat = complex(zeros(b.numstates),zeros(b.numstates));

%% Define Ionic Hamiltonian, H0
% This contains the energies of the basis states relative to the ground states absent of any
% Zeeman shifts, but including the Doppler shift. Note that this includes the rotating wave 
% approximation.
H0 = zeromat;
exc = unique(b.exc);
for k = exc
    H0(k,k) = (v-det);
end

%% Define Zeeman Hamiltonian, Hz
% Calculate the Anomalous Zeeman shift for each excited state
gl = @(j,l,s) 1+(j*(j+1)+s*(s+1)-l*(l+1))/(2*j*(j+1));  % lande g-factor
Hz = zeromat;
for k = 1:b.numstates
    g = gl(b.j(k),b.l(k),b.s(k));
    Hz(k,k) = c.uB*B/c.hbar/c.gam422*g*b.m(k);
end

%% Define Dipole Hamiltonian
Hd = zeromat;
for i = 1:length(b.exc)
    exc = b.exc(i);
    gnd = b.gnd(i);            
    Hd(exc,gnd) = -Om(exc,gnd)/2;
    Hd(gnd,exc) = conj(Hd(exc,gnd));
end

%% Define H, the System Hamiltonian
% Note that H is in units of \hbar\gamma
H = H0+Hz+Hd;

end