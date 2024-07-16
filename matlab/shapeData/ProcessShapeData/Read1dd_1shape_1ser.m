function [dataOut, exists, validFile] = Read1dd_1shape_1ser(serNo, shapeStr, ddStr)
if (nargin < 1)
    serNo = 4999;
end
if (nargin < 2)
    shapeStr = '2';
end
if (nargin < 3)
    ddStr = '0.9';
end

serNo_s = num2str(serNo);
fileName = ['../../../../ShapeData/dd', ddStr, '/shape', shapeStr, '/', serNo_s, '.txt'];
exists = zeros(9, 1);
dataOut = cell(9, 1);
if (~isfile(fileName))
    validFile = 0;
    return;
end
validFile = 1;
dat = load(fileName, '-ascii');
[m, n] = size(dat);
for i = 1:m
    llc = dat(i, 1);
    vec = dat(i, 2:n);
    if (llc < -4)
        lci = 9;
    else
        lci = round(1 - (llc + 0.5) / 0.5);
    end
    dataOut{lci} = vec;
    exists(lci) = 1;
end
