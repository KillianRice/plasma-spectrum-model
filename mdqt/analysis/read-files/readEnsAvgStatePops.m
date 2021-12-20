function [out] = readEnsAvgStatePops(directory)
% directory (string): full path to folder with 'energies.csv'

% read in files from 'directory' and attempt to find file 'energies.csv'
f = filesep;
files = dir(directory); % read in files from 'directory'
ind = strcmp({files.name},'ensAvgStatePops.csv'); % identify file index for 'energies.csv'

% ensure file was found within directory
if ~max(ind)
    error('File ''ensAvgStatePops.csv'' not found within given directory.')
end

% read in contents of 'energies.csv'
opts = delimitedTextImportOptions('Delimiter',',');
opts.VariableTypes = 'double';
statePops = readmatrix([files(ind).folder f files(ind).name],opts);

% transfer contents of 'statePops' to output structure 'out'
if ~isempty(statePops)
    out = struct;
    out.t = statePops(:,1)';
    out.pops = statePops(:,2:2:end-1)';
else
    out = struct;
end


end