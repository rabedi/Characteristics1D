function plot_fragmentation_results(runName, fragmentationType, subdomainNum)
if nargin < 1
    runName = 'run_0';
end
if nargin < 2
    fragmentationType = 1; % 0-> 'D'; 1 -> Max_DelU ; 2 -> DelU 
end
if nargin < 3
    subdomainNum = 0;
end
sd = ['sd_', num2str(subdomainNum), '_'];
preft = '_tAll__FragCrn_';

if (fragmentationType == 0)
    ft = '0_D';
elseif (fragmentationType == 1)
    ft = '1_Max_DelU';
elseif (fragmentationType == 2)
    ft = '2_DelU';
end
rootFolder = ['../_PPS2_', runName];
rt = '';
if (strcmp(rootFolder, '') == 0)
    rt = [rootFolder, '/'];
end

prename = [sd, preft, ft, '_StatFragmentation'];

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
