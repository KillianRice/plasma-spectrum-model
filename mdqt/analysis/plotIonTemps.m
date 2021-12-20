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
datadir = [maindir f 'mdqt' f 'example' f 'dih'];
data = struct;
data.qt = loadRunAvgQT(datadir);
data.md = loadRunAvgMD(datadir);



%% Plot Temperature
c = defineConstants();
ecToK = data.md.pms.Ec*2/c.kB; % dimensionless energy to temperature in K conversion factor

fig = figure;
fig.Position = [490   283   506   420];
fig.Color = [1 1 1];

ax = axes();
hold on
vars = {'Ex','Ey','Ez'};
lgdstr = {'T_x','T_y','T_z'};
l = getLineSpecs(length(lgdstr));

for i = 1:length(vars)
    xdata = [data.md.t].*data.md.pms.wPe*sqrt(3)/2/pi;
    ydata = [data.md.(vars{i})].*ecToK;
    lp = l(i);
    plot(xdata,ydata,lp.style,'LineWidth',2,'MarkerSize',4,'Color',lp.col,'MarkerFaceColor',lp.col,'MarkerEdgeColor',lp.col)
end

ax.PlotBoxAspectRatio = [1 1 1];
ax.FontSize = 12;

xlabel('\omega_p_it/2\pi')
ylabel('T_i (K)')

lgd = legend(lgdstr);
lgd.Position = [0.6275    0.2318    0.1552    0.1405];
