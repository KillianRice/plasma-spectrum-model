%% Initiate Program
clc, clearvars -except inp, close all, f = filesep;

% update Matlab search path with function subdirectories
maindir = extractBefore(matlab.desktop.editor.getActiveFilename,'mdqt');
folders = dir([maindir f 'spectrum-model']);
folders(1:2) = [];
for i = 1:length(folders)
    if folders(i).isdir
        addpath([folders(i).folder f folders(i).name])
    end
end

% add MDQT data processing functions
addpath('read-files')

% specify simulation data directory for analysis (i.e., folder that contains different jobs to be
% averaged together)
datadir = [maindir f 'mdqt' f 'example' f 'lif'];
data = struct;
data.qt = loadRunAvgQT(datadir);
data.md = loadRunAvgMD(datadir);

%% Evaluate RE and REK Model Using Same Conditions as MDQT
c = defineConstants();
b = defineBasisStates();
radpat = integrateRadiationPattern();
G = @(v,sig) 1/sqrt(2*pi*sig^2).*exp(-0.5.*(v./sig).^2);

% define parameters for spectrum model functions
t = [0 data.qt.t(end)]./data.qt.pms.unit.t; % time with units \gamma^-1
det = data.qt.pms.det; % LIF-laser detuning from unperturbed resonance with units \gamma
p0 = mean(data.qt.pops(2,1)); % fraction of ions that begin in mj=+1/2 ground state
P = p0 - (1-p0); % electron-spin polarization of ions
pol = data.qt.pms.pol; % LIF-laser polarization
theta = data.qt.pms.theta; % angle (radians) that projection of B onto x-y plane subtends with y axis of lab coordinate system
phi = 0; % angle (radians) between B and x-y plane of lab coordinate system
gamL = 0; % LIF-laser linewidth - should be zero because MDQT does not account for that
gamD = c.gam1092/c.gam422; % decay rate into D-state
B = data.qt.pms.B; % magnitude of local magnetic field in Tesla
eps = getLaserPolarization(pol,theta,phi,b,B); 
Om0 = data.qt.pms.Om; % Rabi coupling with clebsch gordan coefficient of 1
dEz = getZeemanShift(B,b,c); % Zeeman shift of each basis state
Om = getRabiFrequency(Om0,eps,b); % Rabi coupling included CG coefficient

ecToK = data.md.pms.Ec*2/c.kB; % energy to temperature conversion factor
data.Ti = data.md.Ex(1)*ecToK+.000001; % don't let ion temperature be zero in the model
sigv = sqrt(c.kB*data.Ti/c.mI); % thermal width of velocity distribution in SI
k = 2*pi/c.lam422; % wavenumber for imaging transition
sigd = k*sigv/c.gam422; % Doppler width of velocity distribution with units \gamma
v = linspace(-4*sigd,4*sigd,201); % velocity bins for distribution with units \gamma

% evaluate RE model for each velocity bin and compute ensemble average
data.re = struct;
data.re(length(v)).v = [];
for j = 1:length(v)
    data.re(j).v = v(j);
    data.re(j).G = G(v(j),sigd);
    re = reStatePopModel(t,v(j),det,P,b,Om,dEz,gamL,gamD);
    re.p(5,:) = 1-sum(re.p,1);
    data.re(j).model = re;
    data.re(j).t = data.re(j).model.t;
    data.re(j).pops = re.p;
end
data.reens.t = data.re.t;
data.reens.pops = zeros(size(data.re(1).pops));

vForInt = [data.re.v];
GForInt = [data.re.G];
popsForInt = zeros(length(data.reens.t),length(vForInt));
for k = 1:size(data.re(1).pops,1) % for each state
    for j = 1:length(data.re) % for each velocity
        popsForInt(:,j) = data.re(j).pops(k,:);
    end
    data.reens.pops(k,:) = trapz(vForInt,GForInt.*popsForInt,2);
end
[data.reens.spec,dSdt] = getSpecFromStatePops(data.reens.t,data.reens.pops,radpat,b);
data.reens.dSdt = dSdt;
data.reens.dSdtInt = zeros(size(data.reens.t));
for j = 2:length(data.reens.t)
    data.reens.dSdtInt(j) = trapz(data.reens.t(1:j),dSdt(1:j));
end

% evaluate REK model - inherently computes ensemble average
[data.rek] = rekStatePopModel(t,0,data.Ti,data.md.pms.n,det,P,b,Om,dEz,gamL,gamD,c);
data.rek.p(5,:) = 1 - sum(data.rek.p,1);
[data.rek.spec,dSdt] = getSpecFromStatePops(data.rek.t,data.rek.p,radpat,b);
data.rek.dSdt = dSdt;
data.rek.dSdtInt = zeros(size(data.rek.t));
for j = 2:length(data.rek.t)
    data.rek.dSdtInt(j) = trapz(data.rek.t(1:j),dSdt(1:j));
end

data.n = data.md.pms.n;
data.Gam = getGamma(data.n,data.Ti);

%% Plot Results
fig = figure;
fig.Position = [562    86   665   827];
fig.Color = [1 1 1];

colvar = {''};
rowvar = {''};
row = 3;
col = 2;
num = 5;

ax = cell(row,col);
iter = 0;
lsspec = getLineSpecs(4);
for i = 1:row
    for j = 1:col
        if iter > num - 1, break, end
        iter = iter + 1;
        ax{i,j} = subplot(row,col,iter);
        hold on
        lgdstr = {};

        xdata = [data.qt.t]./data.qt.pms.unit.t;
        ydata = [data.qt.pops(iter,:)];
        plot(xdata,ydata,'.','LineWidth',2,'MarkerSize',15,'Color',lsspec(1).col)
        lgdstr{1} = 'QT';

        xdata = [data.reens.t];
        ydata = [data.reens.pops(iter,:)];
        plot(xdata,ydata,lsspec(2).style,'LineWidth',2,'MarkerSize',15,'Color',lsspec(2).col)
        lgdstr{2} = 'RE';

        xdata = [data.rek.t];
        ydata = [data.rek.p(iter,:)];
        plot(xdata,ydata,lsspec(3).style,'LineWidth',2,'MarkerSize',15,'Color',lsspec(4).col)
        lgdstr{3} = 'REK';

        ax{i,j}.PlotBoxAspectRatio = [1 1 1];
        ax{i,j}.FontSize = 11;

        if i == row, xlabel('\gamma^-^1 t'), end
        ylabel(['|\langle\psi|' num2str(iter) '\rangle|^2'])
    end
end

lgd = legend(lgdstr);
lgd.Position = [0.8325    0.8467    0.1529    0.0733];

ax{i,j} = subplot(row,col,iter+1);
hold on
xdata = [data.qt.t]./data.qt.pms.unit.t;
ydata = [data.qt.dSdt];
plot(xdata,ydata,lsspec(1).style,'LineWidth',2,'MarkerSize',15,'Color',lsspec(1).col)
lgdstr{1} = 'MDQT';

xdata = [data.reens.t];
ydata = [data.reens.dSdt];
plot(xdata,ydata,lsspec(2).style,'LineWidth',2,'MarkerSize',15,'Color',lsspec(2).col)
lgdstr{2} = 'RE';

xdata = [data.rek.t];
ydata = [data.rek.dSdt];
plot(xdata,ydata,lsspec(3).style,'LineWidth',2,'MarkerSize',15,'Color',lsspec(4).col)
lgdstr{3} = 'REK';

ax{i,j}.PlotBoxAspectRatio = [1 1 1];
ax{i,j}.FontSize = 11;

if i == row, xlabel('\gamma^-^1 t'), end
ylabel('dSdt / Density')

an = annotation('textbox');
an.Position = [0.0058    0.8991    0.9919    0.1143];
an.HorizontalAlignment = 'center';
an.VerticalAlignment = 'middle';
an.LineStyle = 'none';
an.FontSize = 12;

dlm = ' - ';
str1 = [data.qt.pms.pol];
str2 = ['theta' num2str(round(data.qt.pms.theta*180/pi,2,'significant'))];
str3 = ['B' num2str(round(data.qt.pms.B*1e4,2,'significant'))];
str4 = ['Om' num2str(round(data.qt.pms.Om,2,'significant'))];
str5 = ['det' num2str(round(data.qt.pms.det,2,'significant'))];
str6 = ['Ti' num2str(round(data.Ti,2,'significant'))];    
str7 = ['n' num2str(round(data.md.pms.n/1e14,2,'significant')) '\times10^1^4'];

an.String = [str1 dlm str2 dlm str3 dlm str4 dlm str5 dlm str6 dlm str7];
