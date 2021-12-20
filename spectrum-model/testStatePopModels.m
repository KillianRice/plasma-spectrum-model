clc, clearvars -except inp, close all, f = filesep;

% update Matlab search path with function subdirectories
maindir = extractBefore(matlab.desktop.editor.getActiveFilename,mfilename());
folders = dir(maindir);
folders(1:2) = [];
for i = 1:length(folders)
    if folders(i).isdir
        addpath([folders(i).folder f folders(i).name])
    end
end

n = 1; % plasma density with units m^-3
B = 10e-4; % magnetic field amplitude in Tesla
v = 0; % hydrodynamic fluid velocity along LIF-laser propagation axis in \gamma/k
tvec = [0 100]; % exposure period with units \gamma^{-1}
P = 0; % electron-spin polarization of the ions [-1 1]
I = 10; % effective LIF-laser intensity in W/m^2
pol = 'lin'; % <left>, <lin>, or <right> - LIF-laser polarization
theta = 40*(pi/180); % angle that projection of B in x-y plane subtends with y axis of lab coordinate system
phi = 0; % angle that B subtends from x-y plane
tol = 1e-6; % tolerance for solving master equations
Ti = .0001; % ion temperature in K
gamL = 0; % LIF-laser linewidth - set to zero here because master equations don't account for it
gamD = 0; % decay rate from 2P_1/2 - 2D_3/2 - set to zero here because master equations don't account for it

[b] = defineBasisStates();
[c] = defineConstants();
[radpat] = integrateRadiationPattern();
[dEz] = getZeemanShift(B,b,c);
[lc] = getLineCenters(v,dEz,b);
[eps] = getLaserPolarization(pol,theta,phi,b,B);
[Om0] = getStateIndepRabiFreq(I,c);
[Om] = getRabiFrequency(Om0,eps,b);

det = lc(1); % detuning of LIF-laser from unperturbed resonance

[fgr] = fgrStatePopModel(tvec,v,det,P,b,Om,dEz,gamL,gamD);
[re] = reStatePopModel(tvec,v,det,P,b,Om,dEz,gamL,gamD);
[rek] = rekStatePopModel(tvec,v,Ti,n,det,P,b,Om,dEz,gamL,gamD,c);
[me] = meStatePopModel(tvec,v,det,B,P,tol,c,b,Om,dEz);

fig = figure;
fig.Position = [389   300   719   523];
fig.Color = [1 1 1];

colvar = {''};
rowvar = {''};
row = 2;
col = 2;
num = 4;

ax = cell(row,col);
iter = 0;
for i = 1:row
    for j = 1:col
        if iter > num - 1, break, end
        iter = iter + 1;
        ax{i,j} = subplot(row,col,iter);
        hold on
        lgdstr = {};

        xdata1 = re.t;
        ydata1 = re.p(iter,:);
        plot(xdata1,ydata1,'-','LineWidth',2,'MarkerSize',15)
        xdata1 = fgr.t;
        ydata1 = fgr.p(iter,:);
        plot(xdata1,ydata1,'-','LineWidth',2,'MarkerSize',15)
        xdata1 = rek.t;
        ydata1 = rek.p(iter,:);
        plot(xdata1,ydata1,'-','LineWidth',2,'MarkerSize',15)
        xdata1 = me.t;
        ydata1 = me.p(iter,:);
        plot(xdata1,ydata1,'-','LineWidth',2,'MarkerSize',15)

        if i == row, xlabel('t (\gamma_4_2_2^-^1)'), end
        ylabel(['p_' num2str(iter) '(t)'])

        ax{i,j}.PlotBoxAspectRatio = [1 1 1];
        ax{i,j}.FontSize = 11;
        ax{i,j}.Position(1) = ax{i,j}.Position(1) - 0;
        ax{i,j}.Position(2) = ax{i,j}.Position(2) - 0;

    end
end
ax{i,j}.PlotBoxAspectRatio = [1 1 1];
ax{i,j}.FontSize = 11;
ax{i,j}.Position(1) = ax{i,j}.Position(1) - 0;
ax{i,j}.Position(2) = ax{i,j}.Position(2) - 0;

lgd = legend({'re','fgr','rek','me'});
lgd.Position = [0.4594    0.4301    0.1001    0.1405];
