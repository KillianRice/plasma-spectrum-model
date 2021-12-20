function [out] = readVelocityDistribution(directory)
% directory (string): full path to folder with 'velBins.csv' and 'velDist.csv'

f = filesep;
files = dir(directory); % read in files from 'directory'

% read 'velDist.csv'
ind = strcmp({files.name},'velDist.csv'); % identify file index for 'velDist.csv'
if ~max(ind), error('File ''velDist.csv'' not found within given directory.'); end
opts = delimitedTextImportOptions('Delimiter',',');
opts.VariableTypes = 'double';
data1 = readmatrix([files(ind).folder f files(ind).name],opts); % SI units

% read 'velBins.csv'
ind = strcmp({files.name},'velBins.csv'); % identify file index for 'velDist.csv'
if ~max(ind), error('File ''velBins.csv'' not found within given directory.'); end
opts = delimitedTextImportOptions('Delimiter',',');
opts.VariableTypes = 'double';
data2 = readmatrix([files(ind).folder f files(ind).name],opts); % SI units

% check to see that contents were read in from file
if ~isempty(data1)
    out.t = data1(:,1); % time in SI units
    out.vel_bins = data2; % velocity bins with SI units
    out.vel_dist = data1(:,2:end); % number of ions per velocity bin
else
    out = struct;
end

end