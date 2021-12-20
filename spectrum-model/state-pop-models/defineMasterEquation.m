function [dpdt] = defineMasterEquation(~,pin,H,gnd,exc,gam,zeromat)
% t (1x1 double): time in units of \gamma^-1 at which rate equations are defined for, not used
% p (32x1 double): contains the real/imag components for each element of the density matrix at time t
% H (4x4 complex double): system Hamiltonian for LIF configuration
% gnd (1x4 double): indices for ground states of decay paths
% gnd (1x4 double): indices for excited states of decay paths
% gam (1x4 double): decay rates (units \gamma) for each decay path
% zeromat (4x4 complex double): complex matrix of zeros, passed into function for speed

%% Reformat p into matrix form
% the first 16 elements of p are real, the second 16 elements of p are real (but indicate the imaginary
% components). 
pvec = complex(pin(1:16),pin(17:32));
% now p(i) contains the ith complex density matrix element, the next few lines format it properly into a
% matrix
numStates = length(gnd);
p = reshape(pvec,numStates,numStates);

%% Define Hgam - decay matrix in master equation
% This section defines the \gamma dependent part of the master equation
H1 = zeromat;
H2 = zeromat;
H3 = zeromat;

for k = 1:length(exc)
    H1(exc(k),:) = H1(exc(k),:) + (gam(k)/2).*p(exc(k),:);
    H2(:,exc(k)) = H2(:,exc(k)) + (gam(k)/2).*p(:,exc(k));
    H3(gnd(k),gnd(k)) = H3(gnd(k),gnd(k)) - (gam(k)/2)*2*p(exc(k),exc(k));
end
Hgam = H1 + H2 + H3;

%% Define Master Equations
dpdt = -1i*(H*p-p*H)-Hgam;

%% Reformat dpdt for output to ode solver
dpdt = dpdt(:);

dpdtR = real(dpdt);
dpdtI = imag(dpdt);

dpdt = [dpdtR;dpdtI];

end