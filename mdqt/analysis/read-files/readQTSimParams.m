function [out] = readQTSimParams(directory)
% directory: (string) path to simulation data folder with file 'qtParams.csv'

% read in files from 'directory' and load information from file named 
% 'qtParamout.csv
f = filesep;
files = dir(directory);
files(1:2) = [];
ind = strcmp({files.name},'qtParams.csv');
opts = delimitedTextImportOptions('Delimiter',',');
opts.VariableTypes = 'string';
qtParams = readmatrix([files(ind).folder f files(ind).name],opts);

out.numStates = str2double(qtParams(1,2));   % number of magnetic sub-levels within simulation
out.mI = str2double(qtParams(2,2));          % mass of Sr+ in kg
out.unit.t = str2double(qtParams(3,2));      % time unit is \gamma^{-1}, all unit.___ in SI units
out.unit.f = str2double(qtParams(4,2));      % frequency unit is \gamma
out.unit.E = str2double(qtParams(5,2));      % energy unit is \hbar\gamma
out.unit.x = str2double(qtParams(6,2));      % length unit is k = 2\pi/\lambda
out.unit.v = str2double(qtParams(7,2));      % velocity unit is (length unit/time unit)
out.N = str2double(qtParams(8,2));           % number of ions in simulation
out.Om = str2double(qtParams(9,2));          % Rabi frequency of LIF imaging laser with units \gamma
out.theta = str2double(qtParams(10,2));      % Angle between imaging polarization and local field in units radians
out.B = str2double(qtParams(11,2));          % magnetic field amplitude in Tesla
out.det = str2double(qtParams(12,2));        % detuning of imaging laser with units \gamma
out.dt = str2double(qtParams(13,2));         % QT time step with units \gamma^{-1}
out.pol = qtParams{14,2};                    % imaging polarization
out.P = str2double(qtParams(15,2));          % initial population of state |2>

end
