function plot_space_spacetime_results(rootFolder, maxTimeStep, step, subdomainNo)
if nargin < 1
    rootFolder = '../2022_12_24_good_fine';
    rootFolder = '../run_0';
end
% if < 0, does not print space files
if nargin < 2
    maxTimeStep = 100000000;
end
if nargin < 3
    step = -1;
end
if nargin < 4
    subdomainNo = 0;
    subdomainNo = 1;
%    subdomainNo = 2;
end
sd = ['sd_', num2str(subdomainNo), '_'];

rt = '';
if (strcmp(rootFolder, '') == 0)
    rt = [rootFolder, '/'];
end


prename = [sd, '_Summary'];
fn = [rt , prename, '.txt'];
fid = fopen(fn, 'r');
if (fid > 0)
    dwh = datawheader;
    dwh.rootFolder = rootFolder;
    dwh = dwh.read(fid);
    fclose(fid);
    dwh.plotRaw(2, prename);
    fclose('all');
    close('all');
end


if (maxTimeStep > 0)
    if (step < 0)
        fid = -1;
        for cntr = 1:10000000
%        for cntr = 1:640000
            fn = [rt , sd, 'BulkInterface_tPos_', num2str(cntr), '.txt'];
            fid = fopen(fn, 'r');
            if fid > 0
                fclose(fid);
                break;
            end
        end
        if (fid > 0)
            step = cntr;
        end
    end
    if (step > 0)
        for cntr = 0:step:maxTimeStep
            prename = [sd, 'BulkInterface_tPos_', num2str(cntr)];
            fn = [rt , prename, '.txt'];
            fid = fopen(fn, 'r');
            if fid > 0
                dwh = datawheader;
                dwh.rootFolder = rootFolder;
                dwh = dwh.read(fid);
                fclose(fid);
            
                dwh.plotRaw(3, prename);
                fclose('all');
                close('all');
            else
                break;
            end
        end
    end
end

