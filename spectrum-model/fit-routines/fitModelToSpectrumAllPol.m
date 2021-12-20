function [out] = fitModelToSpectrumAllPol(dets,spec,specse,x,y,pol,Ti,P,tE,I,model,mag,gamL,weightFit)
% dets (3x1 cell): dets{i} is a laser detuning vector corresponding to pol{i}, units in MHz
% spec (3x1 cell): spec{i} is the spectrum for the corresponding detuning vector, units LIF sig. per unit mm 
% specse (3x1 cell): specse{i} is the spectrum standard error, units LIF sig. per unit mm
% x (double): x-position (mm) of local region
% y (double): y-position (mm) of local region
% pol (3x1 cell): LIF-laser polarization in lab coordinate system ('pi', 'left', or 'right')
% Ti (double): initial guess for local ion temperature in K
% tE (3x1 double): exposure time used to collect spectrum with units \gamma^-1
% I (3x1 double): uniform laser intensity with units W/m^2
% mod (string): select model to fit spectrum with ('fgr','re','me')
% gamL (1x1 double): laser linewidth in MHz
% weightFit (boolean): (true) weight fits by spectrum standard error (false) do not

c = defineConstants(); % load useful constants in SI units
b = defineBasisStates(); % get information for quantum states involved in LIF

% radpat contains the fraction of spontaneously emitted photons we capture using a 1:1 optical
% relay. The captured fraction depends on the angular momentum of the emitted photon.
radpat = integrateRadiationPattern();

%% Normalize Spectrum for Fit
% get normalization factor for spectra
spec_max = zeros(size(spec));
spec_min = zeros(size(spec));
for i = 1:length(pol)
    spec_max(i) = max(spec{i});
    spec_min(i) = min(spec{i});
end
fac = max(spec_max) - min(spec_min);

% normalize spectra to have numbers of order unity
specForFit = cell(size(spec));
specseForFit = cell(size(spec));
for i = 1:length(pol)
    specForFit{i} = spec{i}./fac;
    specseForFit{i} = specse{i}./fac;
end

% use the spectrum's standard error to compute weights for the fit
w = cell(size(specse));
for i = 1:length(pol)
    if weightFit
        w{i} = sqrt(1./specse{i});
        w{i} = w{i}./max(w{i});
    else
        w{i} = ones(size(specse{i}));
    end
end

% convert detunings to units of \gamma, also log number of detunings
detsdim = cell(size(pol));
det_length = zeros(size(pol));
for i = 1:length(pol)
    detsdim{i} = dets{i}.*2*pi*1e6/c.gam422;
    det_length = length(detsdim{i});
end
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

r = 4;
s = cell(size(pol));
specfilt = cell(size(pol));
for i = 1:length(pol)
    % Evaluate spectrum model (with zero velocity)
    s{i} = spectrumModel(detsdim{i},0,g.Ti,g.P,I(i),pol{i},B,theta,phi,tE(i),radpat,model,gamLdim,gamDdim,c);

    specfilt{i} = movingWeightedAverage(specForFit{i},w{i},2,2); % smooth out spectrum for initial guess determination
    spec1 = specfilt{i}-min(specfilt{i}); % ensure spectrum is positive definite
    spec1 = abs(spec1.^r); % create weights for detunings
    spec1 = spec1/sum(spec1);
    det1 = sum(spec1.*detsdim{i}); % compute weighted average detuning

    % compute weighted average of spectrum from model with zero velocity - same process as before
    spec2 = ([s{i}.S]-min([s{i}.S]));
    spec2 = abs(spec2.^r);
    spec2 = spec2/sum(spec2);
    det2 = sum(spec2.*detsdim{i});

    % use the difference in the weighted average detuning as the initial velocity guess
    k = 2*pi/c.lam422; % wavenumber for 422 transition
    g.v(i) = (det2-det1)*c.gam422/k; % velocity guessin units of m/s

    % compute guesses for amplitude and offset
    ind = (specfilt{i}-min(specfilt{i})) < .1.*(max(specfilt{i}) - min(specfilt{i}));
    g.offset(i) = mean(specfilt{i}(ind));
    g.amp(i) = (max(specfilt{i}) - g.offset(i))./max([s{i}.S]);
end

% integrate spectra and use as weights for determining initial guesses
specint = zeros(1,length(pol));
for i = 1:length(pol)
    specint(i) = abs(trapz(detsdim{i},specfilt{i}));
end
specint = specint./sum(specint);
[specint,sortind] = sort(specint);
g.v = g.v(sortind);
g.amp = g.amp(sortind);

% remove element with smallest integrated signal because it is likely the least accurate for a guess
specint(1) = [];
g.v(1) = [];
g.amp(1) = [];

% renormalize spectrum integral weights to sum to one
specint = specint./sum(specint);

g.v = sum(specint.*g.v);
g.amp = abs(sum(specint.*g.amp));


%% Run RE Model Fit
% format data for fit model - each quantity inputted to the model must be a vector
ind = zeros(1,sum(det_length));
xdata = zeros(1,sum(det_length));
ydata = zeros(1,sum(det_length));
wdata = zeros(1,sum(det_length));
iter = 0;
for i = 1:length(pol)
    for j = 1:length(detsdim{i})
        iter = iter + 1;
        ind(iter) = i;
        xdata(iter) = detsdim{i}(j);
        ydata(iter) = specForFit{i}(j);
        wdata(iter) = w{i}(j);
    end
end

% anonymous function for fit model
fun = @(x,xdata) spectrumFitModelAllPol(xdata,ind,x(3),x(4),x(5),x(6),x(2),x(2),x(2),x(1),x(7),I,pol,B,theta,phi,tE,radpat,model,gamLdim,gamDdim,c);
fitmodel = @(x,xdata) fun(x,xdata).*wdata;

% initial Guesses
pms0(1) = g.Ti; lb(1) = g.TiMin; ub(1) = g.TiMax;
pms0(2) = g.v; lb(2) = pms0(2)-75; ub(2) = pms0(2)+75;
pms0(3) = g.amp; lb(3) = pms0(3)/25; ub(3) = pms0(3)*25;
pms0(4) = g.offset(1); lb(4) = pms0(4)-g.amp; ub(4) = pms0(4)+g.amp;
pms0(5) = g.offset(2); lb(5) = pms0(5)-g.amp; ub(5) = pms0(5)+g.amp;
pms0(6) = g.offset(3); lb(6) = pms0(6)-g.amp; ub(6) = pms0(6)+g.amp;
pms0(7) = g.P; lb(7) = -0.25; ub(7) = 1;

% Do fit
options = optimoptions('lsqcurvefit','Algorithm','trust-region-reflective');
options.MaxIterations = 400;

% do fit
[pms,~,R,ef,~,~,J] = lsqcurvefit(fitmodel,pms0,xdata,ydata.*wdata,lb,ub,options);
spec1 = fun(pms,xdata);

% get fit uncertainties
ydataErr = nlpredci(fun,xdata,pms,R,'Jacobian',full(J));    % 95% confidence interval half-widths
ydataErr = abs(ydataErr)*2/3.92;                    % convert half-widths to standard error
ci = nlparci(pms,R,'jacobian',full(J));    % 95% confidence intervals from fit
se = (ci(:,2) - ci(:,1))/3.92;    % convert 95% confidence intervals to standard error

%% Check results of fit

% fig = figure;
% fig.Position = [449         342        1098         361];
% spec2 = fun(pms0,xdata);
% for i = 1:length(pol)
%     ax{i} = subplot(1,3,i);
%     xplot = detsdim{i};
%     
%     plot(xplot,specForFit{i},'.-','LineWidth',1.5,'MarkerSize',16)
%     hold on
%     plot(xplot,spec2(i == ind),'.-','LineWidth',1.5,'MarkerSize',16)
%     plot(xplot,specfilt{i},'.-','LineWidth',1.5,'MarkerSize',16)
%     axes(ax{i})
%     xplot = detsdim{i};
%     hold on
%     plot(xplot,spec1(i == ind),'.-','LineWidth',1.5,'MarkerSize',16)
%     title('RE Fits to Spectrum')
%     xlabel('\Delta (MHz)')  
%     ylabel('LIF Spectrum (a.u.)')
% end
% legend({'data','guess','filt','fit'});
% close(fig)


%% Output Fit Results
dEz = getZeemanShift(B,b,c); % Zeeman shifts of basis states
linec = -(pms(2)*k/c.gam422-(dEz(b.exc)-dEz(b.gnd))); % Zeeman- and Doppler-shifted transition resonance frequencies

for i = 1:length(pol)
    [~,epsM,eps0,epsP] = getLaserPolarization(pol{i},theta,phi,b,B); % spherical tensor components of LIF-laser polarization
    
    out(i).pol = pol{i}; % LIF-laser polarization
    
    out(i).dets = dets{i}; % detuning vector in MHz
    out(i).spec = spec1(i == ind).*fac; % fitted spectrum rescaled to same units as input spectrum
    out(i).specErr = ydataErr(i == ind).*fac; % standard error in spectrum fit rescaled to same units as input spectrum
    out(i).ef = ef; % exit flag from fit

    out(i).Ti = pms(1); % fitted temperature in K
    out(i).v = pms(2); % fitted velocity in m/s
    out(i).amp = pms(3).*fac; % fitted amplitude rescaled to same units as input spectrum
    out(i).offset = pms(3+i).*fac; % fitted offset rescaled to same units as input spectrum
    out(i).TiErr = se(1); % standard error in temperature fit in K
    out(i).vErr = se(2); % standard error of velocity fit in K
    out(i).ampErr = se(3).*fac; % standard error in fit amplitude
    out(i).offsetErr = se(3+i).*fac; % standard error in offset
    out(i).P = pms(7); % electron-spin polarization fit value
    out(i).PErr = se(7); % standard error in P fit value
    
    out(i).I = I(i); % effective local LIF-laser intensity with units W/m^2
    out(i).theta = theta;  % angle of local field relative to y axis in radians
    out(i).phi = phi; % angle that projection of B into x-y plane subtends from y axis in radians
    out(i).B = B; % amplitude of local magnetic field in Tesla
    
    out(i).eps = [epsM eps0 epsP]; % LIF-laser polarization in lab coordinate system
    out(i).dEz = dEz.*c.gam422/2/pi/1e6; % Zeeman shifts of each basis state with units MHz
    out(i).linec = linec.*c.gam422/2/pi/1e6; % Zeeman- and Doppler-shifted transition resonance frequencies with units MHz
end

end


