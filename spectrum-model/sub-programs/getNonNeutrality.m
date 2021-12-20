function [frac] = getNonNeutrality(n0,Te,sig)
%% Notes
% n0 (1x1 double): peak plasma density for exponential plasma (m^-3)
% Te (1x1 double): initial electron temperature (K)
% sig (1x1 double): geometric mean of rms plasma size (m)

% This function calculates the percentage of non-neutrality, expressed as a fraction 'frac' for a
% plasma that is initially exponentially distributed. The x axis is the symmetry axis of the plasma,
% which is the tight axis of the quadrupole magnetic trap. See Thomas Langin PhD thesis page 16 for
% more information.

%% Calculate total number of ions in exponentially distributed plasma
c = defineConstants();

spacelim = sig*100;

x = linspace(-spacelim,spacelim,5001);
y = linspace(-spacelim,spacelim,5001);
z = linspace(-spacelim,spacelim,5001);
n = @(x,y,z) n0.*exp(-sqrt(x.^2+y.^2/4+z.^2/4)./sig);

Ni = integral3(n,min(x),max(x),min(y),max(y),min(z),max(z));

%% Calculate non-neutrality of plasma

% calculate N*
Ns = 3/2*sqrt(pi/2)*(4*pi*c.eps*sig/c.e^2)*c.kB*Te;

% calculate ratio of Ne/Ni = frac
frac = (sqrt(Ni/Ns)-1)/sqrt(Ni/Ns);

end