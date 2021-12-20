function [w_p] = getPlasmaFreq(n,m)
%% Function Parameters

% n (nxm double): plasma density (m^-3)
% m (1x1 double): mass (kg) of species of interest (i.e., electrons or ions)
% w_p (nxm double): plasma oscillation frequency (rad/s)

%% Notes

% This function calculates the Coulomb coupling parameter for a plasma species with density n and
% temperature T

%% Calculate Coulomb Coupling Parameter

c = defineConstants(); % get list of fundamental constants
w_p = sqrt(n.*c.e^2/(m*c.eps)); % plasma frequency (rad/s)

end