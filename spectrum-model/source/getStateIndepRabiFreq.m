function [Om0] = getStateIndepRabiFreq(I,c)
% I (double): LIF laser intensity (W/m^2)
% c (struct): structure with useful constants in SI units
% Om0 (double): sublevel-independent Rabi freq (units gam422)

% This function computes the Rabi coupling induced by the LIF laser on the LIF manifold. Any
% quantities dependent on the magnetic sublevels of a particular transition are included in the full
% Rabi coupling calculation in 'getRabiFrequency.m'.

w = c.c*(2*pi/c.lam422); % angular frequency for the LIF transition
Om0 = sqrt(6*pi*c.gam422*c.c^2*I/(c.hbar*w^3))/c.gam422;

end