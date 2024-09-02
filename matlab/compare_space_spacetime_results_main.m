function compare_space_spacetime_results_main(serialNums, versionNums, subdomainNo, compareTimeVals)
if nargin < 1
    serialNums = [0, 1];
%    serialNums = 0:5;
%    serialNums = [0, 1, 2];
    serialNums = 0;
end
if nargin < 2
%    versionNums = [11, 12, 13];
%    versionNums = 0:1;
%    versionNums = 0;
%    versionNums = 30000:30003;
    versionNums = 40433;
    versionNums = 0;
end
if nargin < 3
    subdomainNo = 0;
%    subdomainNo = 1;
end
if nargin < 4
    compareTimeVals = 1;
end
timeSerialNumbers = [];
timeStamp = '';
compare_space_spacetime_results(1, serialNums, versionNums, timeSerialNumbers, timeStamp, subdomainNo);
if (compareTimeVals == 0)
    return;
end
sd = ['sd_', num2str(subdomainNo), '_'];


sz_serialNums = length(serialNums);
sz_versionNums = length(versionNums);
totalTimeSteps = zeros(sz_serialNums, sz_versionNums);
numTimeStep_BulkInterfacePoints_Print_4PP = zeros(sz_serialNums, sz_versionNums);
maxTime = 0;

cntr = 0;
for sii = 1:sz_serialNums
    si = serialNums(sii);
    sis = num2str(si);
    for vii = 1:sz_versionNums
        vi = versionNums(vii);
        if (vi >= 0)
            vis = num2str(vi);
            vis = ['_V_', vis];
        else
            vis = '';
        end
        cntr = cntr + 1;
        rt = ['../run', vis, '_', sis, '/'];
        prename = ['_', sd, 'keyParameters'];
        fn = [rt , prename, '.txt'];
        fid = fopen(fn, 'r');
        if (fid < 0)
            fprintf(1, 'cannot open file\t%s\n', fn);
            pause
%            exit(0);
        end
        buf = '';
        if (cntr == 1)
            while (strcmp(buf, 'maxTime') == 0)
                buf = fscanf(fid, '%s', 1);
            end
            maxTime = fscanf(fid, '%d', 1);
        end
        while (strcmp(buf, 'totalTimeSteps') == 0)
            buf = fscanf(fid, '%s', 1);
        end
        totalTimeSteps(sii, vii) = fscanf(fid, '%d', 1);
        while (strcmp(buf, 'numTimeStep_BulkInterfacePoints_Print_4PP') == 0)
            buf = fscanf(fid, '%s', 1);
        end
        numTimeStep_BulkInterfacePoints_Print_4PP(sii, vii) = fscanf(fid, '%d', 1);
    end
end

numStp = floor(totalTimeSteps(1, 1) / numTimeStep_BulkInterfacePoints_Print_4PP(1, 1));
for sii = 1:sz_serialNums
    for vii = 1:sz_versionNums
        numStp = min(numStp, floor(totalTimeSteps(sii, vii) / numTimeStep_BulkInterfacePoints_Print_4PP(sii, vii)));
    end
end

for ti = 1:numStp
    timeSerialNumbers = zeros(sz_serialNums, sz_versionNums);
    timeStamp = maxTime / numStp * ti;
    for sii = 1:sz_serialNums
        for vii = 1:sz_versionNums
            timeSerialNumbers(sii, vii) = ti * numTimeStep_BulkInterfacePoints_Print_4PP(sii, vii);
        end
    end
    compare_space_spacetime_results(0, serialNums, versionNums, timeSerialNumbers, timeStamp, subdomainNo);
    fclose('all');
    close('all');
end
