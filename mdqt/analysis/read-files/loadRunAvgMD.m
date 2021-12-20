function [out] = loadRunAvgMD(directory)
% directory (string): full path to simulation data folder - the data folder contains multiple 'run' folders, each of which
% contains .csv files, and are to be averaged together

% read in contents of <directory> - only keep 'run' folder names
f = filesep;
runfolders = dir(directory);
ind = contains({runfolders.name},'run');
runfolders = runfolders(ind);

% load information in each run folder
data = struct;
data(length(runfolders)).run = [];
for i = 1:length(runfolders)
    data(i).run = str2double(extractAfter(runfolders(i).name,'run'));
    data(i).dir = [runfolders(i).folder f runfolders(i).name];
    data(i).opts = readSimOpts(data(i).dir);
    data(i).pms = readPlasmaParams(data(i).dir);
    
    % load kinetic and potential energy data
    E = readEnergyData(data(i).dir);
    varnames1 = fieldnames(E)';
    for j = 1:length(varnames1)
        data(i).(varnames1{j}) = E.(varnames1{j});
    end
    
    % load velocity distribution
    vel = readVelocityDistribution(data(i).dir);
    varnames2 = fieldnames(vel)';
    if ~isempty(varnames2) 
        for j = 1:length(varnames2)
            data(i).(varnames2{j}) = vel.(varnames2{j});
        end
    end
end

% output data averaged over each run
out = struct;
out.opts = data(1).opts;
out.pms = data(1).pms;
variables = unique([varnames1 varnames2]);
for i = 1:length(variables)
    out.(variables{i}) = zeros(size(data(1).(variables{i})));
    for j = 1:length(data)
        out.(variables{i}) = out.(variables{i}) + data(j).(variables{i})./length(data);
    end
end

end