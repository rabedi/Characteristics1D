function [dat, names, num_pts, la, lap] = Periodic_Segment_ReadOneSolve(fileName)

if nargin < 1
    fileName = '../_Solver1D_pfragment__SlnCntr_2_AddedPts0_solveNum_0.txt';
end
fid = fopen(fileName, 'r');
num_c = 117;
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
lap = [];
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
    la(num_pts) = str2num(str{22});
    lap(num_pts) = str2num(str{41});
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
