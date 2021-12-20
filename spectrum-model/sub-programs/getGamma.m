function [Gam] = getGamma(n,T)
%% Function Parameters

% n (nxm double): plasma density with SI units
% T (nxm double): temperature of the species the user wants Gam calculated for (electrons or ions)
%                   in SI units
% Gam (nxm double): Coulomb coupling parameter for plasma species of density n and temperature T

%% Notes

% This function calculates the Coulomb coupling parameter for a plasma species with density n and
% temperature T

%% Calculate Coulomb Coupling Parameter

a = getWignerSeitzRadius(n); % Wigner-Seitz radius in SI units
c = defineConstants(); % get list of fundamental constants
Gam = c.e^2/(4*pi*c.eps*c.kB.*T.*a); % calculate Coulomb coupling parameter

end