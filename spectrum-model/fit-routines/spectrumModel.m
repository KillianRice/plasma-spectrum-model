function [out] = spectrumModel(dets,v,Ti,P,I,pol,B,theta,phi,tE,radpat,model,gamL,gamD,c)
% dets (1xn double): detunings (units gam422) to return spectrum for
% v (double): plasma hydrodynamic flow velocity (m/s)
% Ti (double): ion temperature (K)
% P (double): electron-spin polarization of ions in 2S_{1/2} ground state prior to LIF
% I (double): effective LIF-laser intensity (W/m^2)
% pol (string): LIF-laser polarization in lab coordinate system - 'lin', 'left', or 'right'
% B (double): magnitude of local magnetic field (Tesla)
% theta (double): angle (radians) between y axis of lab coordinate system and projection of local field onto x-y plane
% phi (double): angle (radians) between local magnetic field and x-y plane of lab coordinate system
% tE (double): exposure time (units gam422^-1)
% radpat (struct): fraction of photons captured by optical relay system - see 'integrateRadiationPattern.m'
% model (string): user chooses which model to use ('re' or 'fgr')
% gamL (double): Lorentzian linewidth of LIF laser with units gam422

%% Error Checking
% check that 'model' is specified correctly
if max(strcmp(model,{'re','fgr'})) ~= 1
    error('<model> should be specified as ''re'' or ''fgr''.')
end

%% Initialize Quantum State Information
[t] = [0 tE]; % initial and final time of LIF excitation window (gam422^-1)
[vdim] = v*2*pi/c.gam422/c.lam422;   % convert hydrodynamic flow velocity from m/s to dimensionless units (gam422/k)
[b] = defineBasisStates(); % basis state information for LIF manifold
[eps] = getLaserPolarization(pol,theta,phi,b,B); % LIF laser polarization in local field frame expressed in spherical tensor coordinates
[Om0] = getStateIndepRabiFreq(I,c); % Om0 is Rabi coupling of LIF manifold (independent of sublevel quantum numbers)
[Om] = getRabiFrequency(Om0,eps,b); % Om(i,j) is Rabi coupling between states |i> and |j> induced by LIF laser with units gam422
[dEz] = getZeemanShift(B,b,c); % dEz(k) is the Zeeman shift of transition |k> with units \hbar*gam422
[lc] = getLineCenters(vdim,dEz,b); % linec(t) is the Zeeman- and Doppler- shifted resonance condition each transition

%% Define Detuning Vector to Solve for Spectrum
% This section defines the detuning vector, 'det', for use with the spectrum model. Note the
% difference between 'det', which is defined here, and 'dets', which is the input detuning vector.
% In order to obtain the spectrum for each value in 'dets' we must take a convolution of the LIF
% spectrum lineshape with the local velocity distribution. The input 'dets' does not necessarily 
% have sufficient resolution for the convolution calculation, so here we define the higher resolution
% 'det' vector.

% 'detbgd' holds enough detunings to effectively get the signal baseline
% detline(t,:) holds a more concentrated set of detunings around each transition resonance in 'lc'
% det is the full detuning vector for use with the spectrum model
detbgd = min(dets):max(dets); % one point per 422 transition linewidth
numdets = 21; % number of detunings to use per transition
detline = zeros(numdets,b.numstates);
for i = 1:length(lc)
    detline(:,i) = linspace(lc(i)-1.5,lc(i)+1.5,numdets)';
end
det = unique([detbgd detline(:)']);

%% Evaluate State Pop Model for Each Detuning
fields = {'det','t','p'};
s = cell2struct(cell(length(det),length(fields))',fields);
s(length(det)).det = [];
for j = 1:length(det)
    if strcmp(model,'re')
        s(j) = reStatePopModel(t,vdim,det(j),P,b,Om,dEz,gamL,gamD);
    elseif strcmp(model,'fgr')
        s(j) = fgrStatePopModel(t,vdim,det(j),P,b,Om,dEz,gamL,gamD);
    end
end

%% Calculate LIF Signal from Excited State Populations
for i = 1:length([s.det])
    s(i).S = getSpecFromStatePops(s(i).t,s(i).p,radpat,b);
end

%% Define Maxwell-Boltzmann Velocity Distribution
% The previous section computed the LIF lineshape as a function of laser detuning from unperturbed
% resonance. That treatment of the LIF lineshape is inherently single-particle, meaning the equations 
% are solved for a single ion velocity that is constant in time, and the velocity distribution has not 
% yet been accounted for.

% Assuming local thermal equilibrium, the velocity distribution of the ions is Maxwellian with a
% locally defined temperature. This velocity distribution will inhomogeneously broaden the LIF
% lineshape. Thus, the fluorescence spectrum we measure may be obtained via a convolution of the
% single-particle lineshape obtained in the previous section and the velocity distribution. 

% define width of Maxwell-Boltzmann distribution
sigV = sqrt(c.kB*Ti/c.mI); % RMS width of velocity distribution in m/s
k = 2*pi/c.lam422; % wavenumber for LIF transition
sigG = k*sigV/c.gam422;  % RMS width of velocity distribution with units gam422


convpts = 12; % number of points per linewidth
numLinewidths = 4; % number of linewidths to use when defining detuning vector for convolution
sigTot = numLinewidths*sigG;
totalpts = convpts*numLinewidths*2+1; % odd number of conv points that scales with width
detconv = @(shift) shift+linspace(-sigTot,sigTot,totalpts);               % detunings for convolution
% note that the mean ion velocity is already accounted for in the single-particle lineshape. The distribution 'G'
% is only meant to contain the broadening information, and is therefore centered at zero w.r.t. convolution
G = exp(-(detconv(0)).^2./(2*sigG^2))./sqrt(2*pi*sigG^2); % normalized Maxwell-Boltzmann distribution

Gint = trapz(detconv(0),G);

%% Obtain Fluorescence Spectrum via Convolution with Voigt Broadening
sigconv = zeros(size(dets));
for i = 1:length(dets)
    detsForConv = detconv(dets(i));
    sig = interp1([s.det],[s.S],detsForConv,'linear','extrap');
    sigconv(i) = trapz(detsForConv,sig.*G);
end
sigconv = sigconv./Gint;

out.dets = dets;
out.S = sigconv;

end
