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

x = 2; % position on x axis in mm
y = 0; % position on y axis in mm;
n = 1; % plasma density with units m^-3
v = 30; % hydrodynamic fluid velocity along LIF-laser propagation axis in m/s
tE = 100; % exposure period with units \gamma^{-1}
P = .5; % electron-spin polarization of the ions [-1 1]
I = 100; % effective LIF-laser intensity in W/m^2
pol = 'lin'; % <left>, <lin>, or <right> - LIF-laser polarization
tol = 1e-6; % tolerance for solving master equations
Ti = .1; % ion temperature in K
gamL = c.gamL/c.gam422; % LIF-laser linewidth - set to zero here because master equations don't account for it
gamD = c.gam1092/c.gam422; % decay rate from 2P_1/2 - 2D_3/2 - set to zero here because master equations don't account for it
dets = linspace(-300,300,201)*2*pi*1e6/c.gam422; % LIF-laser detuning in \gamma^{-1}
model = 're';

[~,B,theta,phi] = getQuadrupoleField([x y 0]);
[b] = defineBasisStates();
[radpat] = integrateRadiationPattern();
[dEz] = getZeemanShift(B,b,c);
[lc] = getLineCenters(v,dEz,b);
[eps] = getLaserPolarization(pol,theta,phi,b,B);
[Om0] = getStateIndepRabiFreq(I,c);
[Om] = getRabiFrequency(Om0,eps,b);

[s] = spectrumModel(dets,v,Ti,P,I,pol,B,theta,phi,tE,radpat,model,gamL,gamD,c);

fitopts.weightFits = true;
fitopts.model = model;
fitopts.fitspinpol = true;
[fit] = fitModelToSpectrum(s.dets.*c.gam422/2/pi/1e6,s.S,s.S.*0.1,x,y,pol,Ti,P,tE,I,true,gamL*c.gam422/2/pi/1e6,fitopts);

figure
plot(s.dets,s.S)
hold on
plot(s.dets,fit.spec)
legend({'data','fit'})
