function [out] = fitModelToSpectrum(dets,spec,specse,x,y,pol,Ti,P,tE,I,mag,gamL,fitopts)
% dets (1xn double): laser detuning in MHz
% spec (1xn double): spectrum (counts/volume)
% specse (1xn double): spectrum standard error (counts/volume)
% x (double): x pos (mm) of local region in z = 0 plane
% y (double): y pos (mm) of local region in z = 0 plane
% pol (string): LIF-laser polarization ('lin', 'left', or 'right')
% Ti (double): initial guess for local ion temperature in K
% P (double): initial guess for electron-spin polarization
% tE (double): duration of excitation period used to collect spectrum with units \gam422^{-1}
% I (double): effective LIF-laser intensity with units W/m^2
% mag (boolean): (true) external quadrupole magnetic fields (false) no external magnetic fields
% gamL (double): LIF-laser linewidth in MHz
% fitopts.model (string): select model to fit spectrum with ('fgr' or 're')
% fitopts.weightFits (boolean): (true) weight fits by spectrum standard error (false) fits are not weighted
% fitopts.fitspinpol (boolean): (true) fit for electron-spin polarization (false) fix P = 0

c = defineConstants(); % load useful constants in SI units
b = defineBasisStates(); % get information for quantum states involved in LIF

% radpat contains the fraction of spontaneously emitted photons we capture using a 1:1 optical
% relay. The captured fraction depends on the angular momentum of the emitted photon.
radpat = integrateRadiationPattern();

%% Normalize Spectrum for Fit
% normalize spectra to be of order unity to maximize fit speed
fac = max(spec) - min(spec);
specForFit = spec./fac;
specseForFit = specse./fac;

% derive fit weights from standard error of spectrum
if fitopts.weightFits
    w = sqrt(1./specseForFit);
    w = w./mean(w);
else
    w = ones(size(specseForFit));
end

% convert quantities from units of MHz to units of \gam422^{-1}
detsdim = dets.*2*pi*1e6/c.gam422; % detuning vector
gamLdim = 2*pi*gamL*1e6/c.gam422; % LIF-laser linewidth
gamDdim = c.gam1092/c.gam422; % natural linewidth of 2D_3/2 -> 2P_1/2 transition

%% Set Up Inputs for Models + Obtain Initial Guesses for Fit
% 'g' is a struct that contains initial guesses and lower/upper bounds for fit parameters.

g.Ti = Ti; % guess for ion temperature in K
g.TiMin = .05; % we cannot resolve temperatures lower than this without serious effort
g.TiMax = 5; % upper limit on ion temp because fitted temp can get unphysically high for low SNR spectra
g.P = P; % guess for electron-spin polarization

if mag % plasma is magnetized by quadrupole magnetic fields
    [~,B,theta,phi] = getQuadrupoleField([x y 0]);      % quadrupole magnetic field model (Tesla)
else % no external magnetic fields applied to plasma
    B = 0; 
    phi = nan; % phi is irrelevant for unmagnetized plasma because optical field sets ion quantization axis
    theta = nan; % theta is irrelevant for unmagnetized plasma because optical field sets ion quantization axis
end

% evaluate spectrum model at initial guesses (with zero velocity)
s = spectrumModel(detsdim,0,g.Ti,g.P,I,pol,B,theta,phi,tE,radpat,fitopts.model,gamLdim,gamDdim,c);

% compute weighted average of measured spectrum
r = 4;
specfilt = movingWeightedAverage(specForFit,w,2,2); % smooth out spectrum for determining initial guesses
spec1 = specfilt-min(specfilt); % ensure spectrum is positive definite
spec1 = abs(spec1.^r); % create weights for detunings
spec1 = spec1/sum(spec1); % create weights for detunings
det1 = sum(spec1.*detsdim); % compute weighted average detuning - helps determine velocity initial guess

% compute weighted average of re spectrum model with zero velocity - same process as lines above
spec2 = [s.S]-min([s.S]);
spec2 = abs(spec2.^r);
spec2 = spec2/sum(spec2);
det2 = sum(spec2.*detsdim);

% use the difference in the weighted average detuning as the initial velocity guess
k = 2*pi/c.lam422; % wavenumber for 422 transition
g.v = (det2-det1)*c.gam422/k; % velocity guess in units of m/s - 1st order doppler shift is \delta = k*v

% compute guesses for amplitude and offset
ind = (specfilt-min(specfilt)) < .1.*(max(specfilt) - min(specfilt)); % identify detunings with low signal
g.offset = mean(specfilt(ind)); % offset guess taken as mean of low signal points
g.amp = (max(specfilt) - g.offset)./max([s.S])*1.25; % initial amplitude is taken as peak of filtered signal above the offset guess

%% Fit Spectrum with Model
% data for fit
xdata = detsdim;
ydata = specForFit;

% set up fit model with initial guesses and upper bounds
if mag && fitopts.fitspinpol % spin polarization is fit parameter
    % anonymous function for fit model
    fun = @(x,xdata) spectrumFitModel(xdata,abs(x(3)),x(4),x(2),x(1),x(5),I,pol,B,theta,phi,tE,radpat,fitopts.model,gamLdim,gamDdim,c);
    fitmodel = @(x,xdata) (fun(x,xdata) - ydata).*w; % cost function for weighted fitting

    % initial guesses
    pms0(1) = g.Ti; lb(1) = g.TiMin; ub(1) = g.TiMax; % ion temperature
    pms0(2) = g.v; lb(2) = pms0(2)-75; ub(2) = pms0(2)+75; % hydrodynamic flow velocity
    pms0(3) = g.amp; lb(3) = pms0(3)/25; ub(3) = pms0(3)*25; % spectrum amplitude
    pms0(4) = g.offset; lb(4) = pms0(4) - g.amp; ub(4) = pms0(4) + g.amp;  % spectrum offset
    pms0(5) = g.P; lb(5) = -0.25; ub(5) = 1; % electron-spin polarization
else % electron-spin polarization fixed to initial guess
    % anonymous function for fit model
    fun = @(x,xdata) spectrumFitModel(xdata,abs(x(3)),x(4),x(2),x(1),g.P,I,pol,B,theta,phi,tE,radpat,fitopts.model,gamLdim,gamDdim,c);
    fitmodel = @(x,xdata) (fun(x,xdata) - ydata).*w; % cost function for weighted fitting

    % initial guesses
    pms0(1) = g.Ti; lb(1) = g.TiMin; ub(1) = g.TiMax; % ion temperature
    pms0(2) = g.v; lb(2) = pms0(2)-75; ub(2) = pms0(2)+75; % hydrodynamic flow velocity
    pms0(3) = g.amp; lb(3) = pms0(3)/25; ub(3) = pms0(3)*25; % spectrum amplitude
    pms0(4) = g.offset; lb(4) = pms0(4) - g.amp; ub(4) = pms0(4) + g.amp; % spectrum offset
end

% if the initial guess for the ion temperature is not within range, set it to 100 mK
if ~(g.Ti > lb(1) && g.Ti < ub(1))
    pms0(1) = 0.1;
end

% do fit using standard lsqcurvefit options
options = optimoptions(@lsqcurvefit);
[pms,~,R,ef,~,~,J] = lsqcurvefit(fitmodel,pms0,xdata,zeros(size(ydata)),lb,ub,options);

% get fit uncertainties
ydataErr = nlpredci(fun,xdata,pms,R,'Jacobian',full(J))'; % 95% confidence interval half-widths
ydataErr = abs(ydataErr)*2/3.92; % convert half-widths to standard error
ci = nlparci(pms,R,'jacobian',J); % 95% confidence intervals from fit
se = (ci(:,2) - ci(:,1))/3.92; % convert 95% confidence intervals to standard error

%% Check Fit Results

% plot data, initial guess, and fit result - use for diagnosing initial guess and fit errors
% fig = figure;
% fig.Position = [449   283   570   420];
% errorbar(dets,ydata,specseForFit,'.-','LineWidth',1.5,'MarkerSize',16)
% hold on
% plot(dets,specfilt,'LineWidth',1)
% plot(dets,fun(pms0,xdata),'.-','LineWidth',1.5,'MarkerSize',16)
% hold on
% plot(dets,fun(pms,xdata),'.-','LineWidth',1.5,'MarkerSize',16)
% legend({'data','datafilt','guess','fit'});
% title('Fit to Spectrum')
% xlabel('\Delta (MHz)')
% ylabel('LIF Spectrum (a.u.)')
% close(gcf)

%% Output Fit Results
% all output quantities that are proportional to spectrum amplitude must be rescaled back to OG units using 'fac'

[~,epsM,eps0,epsP] = getLaserPolarization(pol,theta,phi,b,B); % laser polarization in local field frame using spherical tensor notation
dEz = getZeemanShift(B,b,c); % zeeman shift for each basis state
linec = getLineCenters(pms(2)*k/c.gam422,dEz,b); % Zeeman- and Doppler-shifted resonance condition for each transition

out.dets = dets; % detuning vector in MHz
out.spec = fun(pms,detsdim).*fac;    % fitted spectrum rescaled to same units as inputted spectrum
out.specErr = ydataErr.*fac;    % standard error in spectrum fit rescaled to same units as inputted spectrum
out.ef = ef; % error flag from fit

out.Ti = pms(1); % fitted temperature in K
out.v = pms(2); % fitted velocity in m/s
out.amp = pms(3)*fac; % fitted amplitude rescaled to same units as input spectrum
out.offset = pms(4)*fac; % fitted spectrum offset rescaled to same units as input spectrum
out.TiErr = se(1);  % standard error in Ti in K
out.vErr = se(2); % standard error in v in m/s
out.ampErr = se(3).*fac; % standard error in spectrum amplitude
out.offsetErr = se(4).*fac; % standard error in spectrum offset
if mag && fitopts.fitspinpol % if spin pol used as fit parameter
    out.P = pms(5); % spin pol resulting from fit
    out.PErr = se(5);
else % spin pol. not used as fit parameter
    out.P = g.P; % spin pol fixed to initial value
    out.PErr = 0; % no error indicates it wasn't fit
end

out.I = I;    % LIF-laser intensity W/m^2
out.theta = theta;  % angle of local field relative to y axis with units radians
out.phi = phi; % angle projection of B into x-y plane subtends from y axis
out.B = B;  % magnetic field in Tesla

out.eps = [epsM eps0 epsP]; % spherical tensor components of laser polarization in the local field frame
out.dEz = dEz.*c.gam422/2/pi/1e6; % Zeeman shifts for each transition with units MHz
out.linec = linec.*c.gam422/2/pi/1e6; % Zeeman/Doppler-shifted resonance condition for each transition in MHz

end
