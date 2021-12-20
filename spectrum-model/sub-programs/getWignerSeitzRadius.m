function [a_ws] = getWignerSeitzRadius(n)
%% Function Parameters

% n (nxm double): plasma density in SI units

% a (nxm double): a(i,j) is the Wigner-Seitz radius (SI units) for a plasma with density n(i,j) (in
% SI units)

%% Notes

% The Wigner-Seitz radius is the average separation of particles in a plasma with density n.

%% Calculate Wigner Seitz Radius for Plasma of Density n

a_ws = (3./(4*pi.*n)).^(1/3); % SI units

end