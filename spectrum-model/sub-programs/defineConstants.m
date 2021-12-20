function [c] = defineConstants()
% This function defines several constants that are generally useful in calculations involving
% Sr+ ultracold neutral plasmas. Some are fundamental constants and others are specific to
% Sr+. All constants are defined in SI units
    
c.mE = 9.1095e-31;              % electron mass in kg
c.kB = 1.3806e-23;              % Boltzmann constant in SI units
c.eps = 8.8542e-12;             % vacuum permittivity in SI units
c.u0 = 1.2566e-6;               % vacuum permeability in SI units
c.uB = 9.274e-24;               % Bohr magneton in J/T
c.mI = 1.455e-25;               % Strontium ion mass in kg
c.h = 6.6261e-34;               % Planck's constant in SI units
c.hbar = c.h/2/pi;              % hbar in SI units
c.c = 299792458;                % speed of light in m/s
c.e = 1.602e-19;                % electron charge in C

c.gam422 = 2*pi*20.3e6;         % natural linewidth of 2S_1/2 -> 2P_1/2 imaging transition
c.gam408 = 2*pi*22.4e6;         % natural linewdith of 2S_1/2 -> 2P_3/2 cooling transition
c.gam1033 = 2*pi*1.384e6;       % natural linewidth of 2D_5/2 -> 2P_3/2 repump transition
c.gam1092 = 2*pi*1.4324e6;      % natural linewidth of 2D_3/2 -> 2P_1/2 repump transition
c.gamL = 2*pi*5e6;              % 422 laser linewidth

c.lam422 = 421.552e-9;          % wavelength of imaging transition in m
c.lam408 = 407.771e-9;          % wavelength of cooling transition in m
c.lam1033 = 1032.7311e-9;       % wavelength of strong repump transition in m
c.lam1092 = 1091.4887e-9;       % wavelength of weak repump transition in m

c.Isat422Pi = 1.0626e3;         % saturation intensity of 422 pi imaging transition in W/m^2
c.Isat422Sig = 531.277;         % saturation intensity of 422 sigma imaging transitions in W/m^2        

end