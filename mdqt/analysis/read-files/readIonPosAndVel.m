function [out] = readIonPosAndVel(filepath)
% directory (string): full path to simulation data folder with file 'ionPosAndVel.csv'

%% Function Notes
% This function loads the user-defined MDQT simulation options. These are things that are not associated with
% either the MD or QT class. The only numerical quantity loaded by this function is the ion temperature, Ti,
% which is in units of K.

%% Function
f = filesep;
files = dir(directory);
ind = strcmp({files.name},'ionPosAndVel.csv');
opts = delimitedTextImportOptions('Delimiter',',');
opts.VariableTypes = 'double';
out = struct;
out.pv = readmatrix([files(ind).folder f files(ind).name],opts);
out.pmean = mean(out.pv(:,1:3),1);
out.vmean = mean(out.pv(:,4:6),1);

end
