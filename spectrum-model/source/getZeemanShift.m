function [dEz] = getZeemanShift(B,b,c)
% B (double): amplitude of local magnetic field (in Tesla)
% b (struct): contains information for basis states involved in LIF, see 'defineBasisStates.m'
% c (struct): contains useful constants for calculations, see 'defineConstants.m'
% dEz (1x4 double): dEz(k) is the Zeeman shift (with units \hbar\gamma) of basis state |k> 

% This function calculates the Anomalous Zeeman shift for each basis state within 'b' for a
% magnetic field of amplitude 'B'.

dEz = zeros(1,b.numstates);
for k = 1:b.numstates
    g = getLandeGFac(b.j(k),b.l(k),b.s(k));
    dEz(k) = g*b.m(k)*c.uB*B/c.hbar/c.gam422;
end

end