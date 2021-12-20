function [lambda] = getDebyeLength(n,Te)
%% Function Parameters

% n (nxm double): plasma density with SI units
% Te (nxm double): electron temperature with SI units
% lambda (nxm double): plasma debye length with SI units

%% Notes

% This function calculates the electron debye length of a neutral plasma with density n and electron
% temperature Te

%% Calculate Coulomb Coupling Parameter

c = defineConstants(); % get list of fundamental constants
lambda = sqrt(c.eps*c.kB.*Te./(n.*c.e^2)); % calculate Coulomb coupling parameter

end