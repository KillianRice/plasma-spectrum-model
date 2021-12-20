function [out] = loadRunAvgQT(directory)
% directory (string): full path to simulation data folder - the data folder contains multiple 'run' folders, each of which
% contains .csv files, and are to be averaged together

[radpat] = integrateRadiationPattern();
b = defineBasisStates();

% read in contents of <directory> - only keep 'run' folder names
f = filesep;
runfolder = dir(directory);
ind = contains({runfolder.name},'run');
runfolder = runfolder(ind);

% load information in each run folder
data = struct;
data(length(runfolder)).run = [];
for i = 1:length(runfolder)
    % record run number and directory
    data(i).run = str2double(extractAfter(runfolder(i).name,'run'));
    data(i).dir = [runfolder(i).folder f runfolder(i).name];
    data(i).basis = readBasisStates(data(i).dir);
    data(i).pms = readQTSimParams(data(i).dir);

    % load state populations
    [sp] = readEnsAvgStatePops(data(i).dir);
    varnames = fieldnames(sp);
    for j = 1:length(varnames)
        data(i).(varnames{j}) = sp.(varnames{j});
    end

    % compute LIF spectrum per unit density - only do if ens avg state pops are recorded
    if ~isempty(varnames)
        [data(i).spec,data(i).dSdt] = getSpecFromStatePops(data(i).t,data(i).pops,radpat,b);
    end
end

% output data averaged over each run
out = struct;
out.basis = data(1).basis;
out.pms = data(1).pms;
if ~isempty(varnames)
    out.t = zeros(size(data(1).t));
    out.pops = zeros(size(data(1).pops));
    out.dSdt = zeros(size(data(1).t));
    out.dSdtInt = zeros(size(data(1).t));
    out.spec = 0;
    for i = 1:length(data) % for each job
        out.t = out.t + data(i).t./length(data);
        out.pops = out.pops + data(i).pops./length(data);
        out.dSdt = out.dSdt + data(i).dSdt./length(data);
        out.spec = out.spec + data(i).spec./length(data);
    end

    out.dSdtInt(1) = 0;
    for i = 2:length(out.t)
        out.dSdtInt(i) = trapz(out.t(1:i)/out.pms.unit.t,out.dSdt(1:i));
    end
end
end
