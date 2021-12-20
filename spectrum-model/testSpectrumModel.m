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
[c] = defineConstants();

n = 1; % plasma density with units m^-3
B = 50e-4; % magnetic field amplitude in Tesla
v = 30; % hydrodynamic fluid velocity along LIF-laser propagation axis in m/s
tE = 255; % exposure period with units \gamma^{-1}
P = .5; % electron-spin polarization of the ions [-1 1]
I = 100; % effective LIF-laser intensity in W/m^2
pol = 'lin'; % <left>, <lin>, or <right> - LIF-laser polarization
theta = 30*(pi/180); % angle that projection of B in x-y plane subtends with y axis of lab coordinate system
phi = 0; % angle that B subtends from x-y plane
tol = 1e-6; % tolerance for solving master equations
Ti = .1; % ion temperature in K
gamL = c.gamL/c.gam422; % LIF-laser linewidth - set to zero here because master equations don't account for it
gamD = c.gam1092/c.gam422; % decay rate from 2P_1/2 - 2D_3/2 - set to zero here because master equations don't account for it
dets = linspace(-300,300,201)*2*pi*1e6/c.gam422; % LIF-laser detuning in \gamma^{-1}

[b] = defineBasisStates();
[radpat] = integrateRadiationPattern();
[dEz] = getZeemanShift(B,b,c);
[lc] = getLineCenters(v,dEz,b);
[eps] = getLaserPolarization(pol,theta,phi,b,B);
[Om0] = getStateIndepRabiFreq(I,c);
[Om] = getRabiFrequency(Om0,eps,b);

[re] = spectrumModel(dets,v,Ti,P,I,pol,B,theta,phi,tE,radpat,'re',gamL,gamD,c);
[fgr] = spectrumModel(dets,v,Ti,P,I,pol,B,theta,phi,tE,radpat,'fgr',gamL,gamD,c);

fig = figure;
fig.Position = [549   315   560   420];
lines = getLineSpecs(2);
plot(dets,fgr.S,lines(1).style,'Color',lines(1).col,'LineWidth',2)
hold on
plot(dets,re.S,lines(2).style,'Color',lines(2).col,'LineWidth',2)

legend({'fgr','re'})