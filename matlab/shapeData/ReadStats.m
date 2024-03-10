function [sNs, mappedStats, snStats] = ReadStats(llc, shape, dd2)
if (nargin < 1)
    llc = -1;
end
if (nargin < 2)
    shape = 1.5;
end
if (nargin < 3)
    dd2 = 0.5;
end
llcs = num2str(llc);
add1 = (llc < -4.0);
shapes = num2str(shape);
dd2s = num2str(dd2);
fn = ['D:\shapeData\llc_',llcs, '/fnb_shape_', shapes, '_llc_', llcs, '_dd2_', dd2s, '_stat.csv'];
fid = fopen(fn, 'r');
if (fid < 0)
    fprintf(1, 'Cannot open file %s\n', fn);
end
tline = fgetl(fid);
sNs = [];
mappedStats = [];
snStats = [];

while (true)
    [buf, neof] = fscanf(fid, '%f', 1);
    if (~neof)
        return;
    end
    buf = fscanf(fid, '%c', 1);
    buf = fscanf(fid, '%f', 1);
    buf = fscanf(fid, '%c', 1);
    buf = fscanf(fid, '%f', 1);
    buf = fscanf(fid, '%c', 1);
    serN = fscanf(fid, '%d', 1);
    row = serN + 1;
    sNs(row) = serN;
    for j = 1:19
        buf = fscanf(fid, '%c', 1);
        tmp(j) = fscanf(fid, '%f', 1);
    end
    if (add1)
        tmp(14) = tmp(14) + 1;
    end
    mappedStats(row, :) = tmp;
    for j = 1:19
        buf = fscanf(fid, '%c', 1);
        tmp(j) = fscanf(fid, '%f', 1);
    end
    if (add1)
        tmp(14) = tmp(14) + 1;
    end
    snStats(row, :) = tmp;
end
fclose(fid);