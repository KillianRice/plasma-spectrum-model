function [out] = readVelocityBins(directory)
% directory (string): full path to folder with 'velBins.csv'

%% Load Ensemble-Averaged Particle velBins
% read in contents of 'velBins.csv'. 
% this file contains the ensemble-averaged kinetic and potential energy as a function of time

% read in files from 'directory' and attempt to find file 'velBins.csv'
f = filesep;
files = dir(directory); % read in files from 'directory'
ind = strcmp({files.name},'velBins.csv'); % identify file index for 'velBins.csv'

% ensure file was found within directory
if ~max(ind)
    error('File ''velBins.csv'' not found within given directory.')
end

% read in contents of 'velBins.csv'
opts = delimitedTextImportOptions('Delimiter',',');
opts.VariableTypes = 'double';
vbins = readmatrix([files(ind).folder f files(ind).name],opts); % SI units

% check to see that contents were read in from file
if isempty(vbins)
    error('No contents were read in from ''velBins.csv''. Maybe the file is empty.')
end

out = vbins;

