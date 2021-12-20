function [out] = readBasisStates(directory)
% directory (string): full path to folder containing 'basisStates.csv'

% read in contents of 'basisStates.csv', which contains the quantum numbers for each quantum state 
% used in the QT simulation.
f = filesep;
files = dir(directory); % read in files from 'directory'
ind = strcmp({files.name},'basisStates.csv'); % identify file index for 'basisStates.csv'

% ensure file was found within directory
if ~max(ind)
    error('File ''basisStates.csv'' not found within given directory.')
end

% read from file
opts = delimitedTextImportOptions('Delimiter',',');
opts.VariableTypes = 'double';
basis = readmatrix([files(ind).folder f files(ind).name],opts);

% check to see that contents were read in from file
if isempty(basis)
    error('No contents were read in from ''basisStates.csv''. Maybe the file is empty.')
end

% reformat matrix 'basis' into output structure 'out'
out = struct; % initialize output structure
out(size(basis,1)).n = []; % initialize size of output structure
variables = {'n','l','s','j','m','manifold'}; % quantum states contained in file
for i = 1:length(out) % for each basis state
    for j = 1:length(variables) % for each quantum number
        out(i).(variables{j}) = basis(i,j);
    end
end

end
