function [kappa] = getScreeningParameter(n,Te)
%% Function Parameters

% n (nxm double): plasma density with SI units
% Te (nxm double): electron temperature with SI units
% kappa (nxm double): plasma screening parameter (dimensionless)

%% Notes

% This function calculates the Coulomb coupling parameter for a plasma species with density n and
% temperature T

%% Calculate Coulomb Coupling Parameter

a = getWignerSeitzRadius(n); % Wigner-Seitz radius in SI units
lambda = getDebyeLength(n,Te); % Debye length with SI units
kappa = a./lambda; % plasma screening parameter

end