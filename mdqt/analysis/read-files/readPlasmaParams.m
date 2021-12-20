function [out] = readPlasmaParams(directory)
% directory (string): full path to folder with 'plasmaParams.csv'

f = filesep;
files = dir(directory);
ind = strcmp({files.name},'plasmaParams.csv');
params = readmatrix([files(ind).folder f files(ind).name]);

% extract each parameter from matrix that was read in from file
out.N = params(1,2);      % number of ions in simulation
out.n = params(2,2);      % constant, uniform plasma density (m^{-3})
out.Ge = params(3,2);     % electron Coulomb coupling parameter
out.L = params(4,2);      % length of simulation box with units a
out.wPe = params(5,2);    % plasma frequency (s^{-1}
out.a = params(6,2);      % Wigner-Seitz radius (m)
out.Ec = params(7,2);     % Ec = e^2/(4*\pi*\eps0*a), natural unit of energy with units J
out.Te = params(8,2);     % electron temperature in K
out.lDeb = params(9,2);   % Debye screening lengh with units a
out.kap = params(10,2);   % plasma screening parameter k = 2*\pi/\lambda (dimensionless)
out.dt = params(11,2);    % md timestep with units w_{pE}^{-1}

end
