function [out] = rekStatePopModel(t,vc,Ti,n,det,P,b,Om,dEz,gamL,gamD,c)
% t (1x2 double): time window to return state pop solutions in units of gamma^-1: t(1) and t(2) are start and end times, respectively
% vc (double): hydrodynamic flow velocity with units gam422/k
% Ti (double): ion temperature in K
% n (double): plasma density in m^-3
% det (double): LIF-laser detuning from unperturbed resonance with units gam422
% P (double): electron-spin polarization
% b (struct): basis state information, see 'basis/defineBasisStates.m'
% eps (double): eps(i,j) is polarization for transition |j> to |i>
% Om (4x4 complex double): Rabi coupling with units gam422
% dEz (1x4 double): dEz(k) is the Zeeman shift (with units \hbargam422) of basis state |k> 
% gamL (double): LIF-laser linewidth in units gam422
% gamD (double): decay rate into 2D_{3/2} state units gam422

[mu] = getCollisionalRelaxationRateMD(n,Ti)/c.gam422; % collisional relaxation rate with units gam422 (see PhysRevLett.110.235001)

%% Define Maxwell-Boltzmann Velocity Distribution
sigV = sqrt(c.kB*Ti/c.mI); % RMS width of velocity distribution in m/s
k = 2*pi/c.lam422; % wavenumber for 422-nm transition 
sigD = k*sigV/c.gam422; % RMS width of velocity distribution in frequency units of \gamma
v = linspace(-4*sigD,4*sigD,101); % vector of velocities spanning entire distribution
G = exp(-(v).^2/(2*sigD^2))/sqrt(2*pi*sigD^2);    % Maxwell-Boltzmann distribution

%% Define Scattering Rate for Velocity Class
% Det(k) is the Doppler- and Zeeman-shifted resonance condition with units \gamma for velocity v(k)
% R(i,j,k) is the scattering rate (with units \gamma) from state |j> to state |i> for velocity v(k)
R = zeros(b.numstates,b.numstates,length(v));
for i = 1:length(v)
    R(:,:,i) = getScatteringRate(det,vc+v(i),dEz,Om,gamL,b);
end

%% Solve REK Equations for Ion State Populations as a Function of Velocity
% gam(i,j) is the decay rate for transition |i> to |j>
gam = zeros(b.numstates);
for k = 1:b.numstates
    gam(b.exc(k),b.gnd(k)) = b.gam(k);
end

% set up column vector for initial state population distribution functions
pInit = [(1-P)/2; (1+P)/2; 0; 0].*G;
pInit = pInit(:);

[w] = getSimpsonWeights(v);
[dpdt] = @(t,p) defineRateEquationsWithKrook(t,p,R,gam,G,w,mu,gamD);

opts = odeset('RelTol',1e-6); % specify solver options

% solve rate equations
[tSol,pSol] = ode45(dpdt,t,pInit,opts);

% reformat
s = struct;
s(length(tSol)).t = [];
for i = 1:length(tSol) % for each time point
    s(i).t = tSol(i);
    s(i).p = reshape(pSol(i,:)',b.numstates,[]);
    s(i).pEns = s(i).p*w;
end

out.det = det;
out.t = tSol';
out.p = [s.pEns];

end