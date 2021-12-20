function [out] = readEnergyData(directory)
% directory (string): full path to folder with 'energies.csv'

f = filesep;
files = dir(directory); % read in files from 'directory'
ind = strcmp({files.name},'energies.csv'); % identify file index for 'energies.csv'
if ~max(ind), error('File ''energies.csv'' not found within given directory.'); end % throw error if not found

opts = delimitedTextImportOptions('Delimiter',',');
opts.VariableTypes = 'double';
E = readmatrix([files(ind).folder f files(ind).name],opts);
if isempty(E)
    error('No contents were read in from ''energies.csv''. Maybe the file is empty.')
end

out.t = E(:,1);       % time t[i] in SI units
out.Ex = E(:,2);      % ensemble-averaged kinetic energy Ex[i] (units Ec) along the x-axis
out.Ey = E(:,3);      % ensemble-averaged kinetic energy Ey[i] (units Ec) along the y-axis
out.Ez = E(:,4);      % ensemble-averaged kinetic energy Ez[i] (units Ec) along the z-axis
out.Ep = E(:,5);      % ensemble-averaged Yukawa potential energy Ep[i] (units Ec)
out.Etot = E(:,6);    % Etot[i] = Ex[i]+Ey[i]+Ez[i]+Ep[i] (units Ec)

end