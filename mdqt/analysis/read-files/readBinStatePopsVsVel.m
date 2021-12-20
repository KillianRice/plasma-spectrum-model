function [out] = readBinStatePopsVsVel(directory)
% directory (string): full path to folder with 'binStatePopsVsVel.csv'

%% Load Ensemble-Averaged Particle velBins
% read in contents of 'binStatePopsVsVel.csv'. 
% this file contains the ensemble-averaged kinetic and potential energy as a function of time

% read in files from 'directory' and attempt to find file 'binStatePopsVsVel.csv'
f = filesep;
files = dir(directory); % read in files from 'directory'
ind = strcmp({files.name},'binStatePopsVsVel.csv'); % identify file index for 'binStatePopsVsVel.csv'

% ensure file was found within directory
if ~max(ind)
    error('File ''binStatePopsVsVel.csv'' not found within given directory.')
end

% read in contents of 'binStatePopsVsVel.csv'
opts = delimitedTextImportOptions('Delimiter',',');
opts.VariableTypes = 'double';
data = readmatrix([files(ind).folder f files(ind).name],opts); % SI units

% check to see that contents were read in from file
if isempty(data)
    error('No contents were read in from ''binStatePopsVsVel.csv''. Maybe the file is empty.')
end

out = struct;
out.t = unique(data(:,1)); % time in SI units
out.numStates = length(find(out.t(1) == data(:,1)));
out.numBins = size(data,1)/out.numStates;
out.pdist = struct;
out.pdist(length(out.t)).t = [];
for i = 1:length(out.t)
    ind = find(data(:,1) == out.t(i));
    for j = 1:length(ind)
        out.pdist(i).t = out.t(i);
        out.pdist(i).(['p' num2str(j)]) = data(ind(j),2:end);
    end
end

end