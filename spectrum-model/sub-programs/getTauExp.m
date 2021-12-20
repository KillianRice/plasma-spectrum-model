function [tau_exp] = getTauExp(sig,Te0,usebeta,beta)
%% Function Parameters

% sig (1x1 double): geometric mean of RMS plasma size in SI units
% Te0 (1x1 double): initial electron temperature with SI units
% tau_exp (1x1 double): hydrodynamic expansion timescale in SI units
% usebeta (bool): (true) use cuspy defintion (false) use gaussiand efinition

if nargin < 3
    usebeta = false;
end

if usebeta && nargin < 4
    beta = 0.63;
else
    beta = 1;
end

%% Notes

% This function calculates the hydrodynamic expansion timescale for a plasma. This function can
% calculate the expansion for either a Gaussian (beta = false) or an Exponential (beta = true)
% plasma. For more information on how the expansion of an exponential plasma differs from a
% Gaussian, see https://arxiv.org/abs/2012.14577.

%% Calculate Coulomb Coupling Parameter
c = defineConstants();



tau_exp = sqrt(c.mI*sig^2/(c.kB*Te0))*beta;