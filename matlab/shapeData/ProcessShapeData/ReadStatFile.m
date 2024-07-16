function [vals, names, n] = ReadStatFile(fileName)
if (nargin < 1)
    fileName = 'stat/fld2/shape1.5_dd0.1_llc-0.5_fld2_stat.txt';
end
fid = fopen(fileName, 'r');
if (fid < 0)
    fprintf(1, 'cannot open file %s\n', fileName);
    pause;
end

% Read the first line from the file
firstLine = fgetl(fid);
% Split the first line into individual strings
names = strsplit(firstLine);
n = length(names);
for i = 1:n
    vals(i) = fscanf(fid, '%g', 1);
end
