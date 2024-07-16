function [dat, names, num_pts, la, la2p, laGp] = Periodic_Segment_ReadOneSolve(fileName)

if nargin < 1
    fileName = '../_Solver1D_pfragment__SlnCntr_2_AddedPts0_solveNum_0.txt';
end
fid = fopen(fileName, 'r');
if (fid < 0)
    fprintf(1, 'fileName %s cannot be opened\n', fileName);
end

num_c = 130; %117;
% 2024/07/15
num_c = 140;
offset = 20;
len = num_c - offset;
for i = 1:num_c
    str{i} = fscanf(fid, '%s', 1);
end
for i = 1:len
    tmp = str{i + offset};
    tmp = strrep(tmp,'/','');
    tmp = strrep(tmp,'\','');
    names{i} = tmp;
end
dat = [];
num_pts = 0;
la = [];
la2p = [];
laGp = [];
str = [];
for i = 1:num_c
    str{i} = fscanf(fid, '%s', 1);
    if (feof(fid))
        break;
    end
end
while (length(str) == num_c)
    num_pts = num_pts + 1;
    valid = strcmp(str{3}, 'yes');
    for i = 1:len
        str_v = str{i + offset};
        [val, tf] = str2num(str_v);
        if (tf == 0)
            val = nan;
        end
        dat(num_pts, i) = val;
    end
    la(num_pts) = str2num(str{23}); %22
    la2p(num_pts) = str2num(str{29}); %41
    laGp(num_pts) = str2num(str{35});
    if (~valid)
        for i = 1:len
            dat(num_pts, i) = nan;
        end
    end
    str = [];
    for i = 1:num_c
        str{i} = fscanf(fid, '%s', 1);
        if (feof(fid))
            break;
        end
    end
end

fclose('all');
close('all');

[la, inds] = sort(la);
b_ordered = 1;
for i = 1:num_pts
    if (inds(i) ~= i)
        b_ordered = 0;
        break;
    end
end
if (b_ordered == 1)
    return;
end
la2pBK = la2p;
laGpBK = laGp;
datBK = dat;
for i = 1:num_pts
    ind = inds(i);
    la2p(i) = la2pBK(ind);
    laGp(i) = laGpBK(ind);
    dat(i, :) = datBK(ind, :);
end
