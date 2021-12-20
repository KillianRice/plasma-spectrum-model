function [out] = readSimOpts(directory)
% directory (string): full path to 'options.csv'

f = filesep;
files = dir(directory);
ind = strcmp({files.name},'options.csv');
options = readcell([files(ind).folder f files(ind).name]);

for rowInd = 1:size(options,1)
    out.(options{rowInd,1}) = options{rowInd,2};
end

end
