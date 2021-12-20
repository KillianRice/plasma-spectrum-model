function [lc] = getLineCenters(v,dEz,b)
% v (double): hydrodynamic flow velocity in gam422/k where k = 2*pi/lam422
% dEz (1x4 double): dEz(k) is the Zeeman shift of state |k> with units \hbar*gam422
% lc (1x4 double): resonance condition for each transition, of which there are 4 driven by LIF laser 

% This function computes the resonance frequency of each Doppler/Zeeman-shifted transition in 
% the 2S_{1/2} - 2P_{1/2} LIF manifold (i.e., the sigma +\- and both pi transitions). We use the
% first order Doppler shift and the anamolous Zeeman shift.

% Note that b.exc(i) and b.gnd(i) contain the excited and ground states for the ith transition,
% respectively. 

lc = dEz(b.exc)-dEz(b.gnd)-v;

end