function TranslateCDFiles(alpha, beta, log2NumSegments, numRealizations)
if nargin < 1
    alpha = '1.95';
end
if nargin < 2
    beta = '0.1';
end
if nargin < 3
    log2NumSegments = 16;
end
if nargin < 4
    numRealizations = 1000;
    numRealizations = 200;
end
nSeg = 2^log2NumSegments;
nSegSer = num2str(nSeg);
nPts = nSeg + 1;
nPtsSer = num2str(nPts);

rootFolder = ['../InhomogeneousFiles/CD/al', alpha, 'b', beta];
targetFolder = [rootFolder, '_np', nPtsSer];
%targetFolder = [rootFolder, '_np', nSegSer];
[status,msg,msgID] = mkdir(targetFolder);
for i = 1:numRealizations
    ii = i - 1;
    istr = num2str(i, '%0.3f');
    fn = [rootFolder, '/cauchyaa', nSegSer, '_al', alpha, '_b', beta, '_', istr, '.csv'];
    fidi = fopen(fn, 'r');
    if (fidi < 0)
        return;
    end
    fclose(fidi);
%    datVals = cell(nPts, 1);
%    for j = 1:nPts
%       datVals{j} = fscanf(fidi, '%f');
%   end
%   fclose(fidi);
    dt = load(fn, '-ascii');
    [nr, nc] = size(dt);
    if nr ~= nSeg
        fprintf(1, 'nr (%d) ~= nSeg (%d)\n', nr, nSeg);
        pause;
    end
    targetName = [targetFolder, '/initial_values_', num2str(ii), '.txt'];
    fid = fopen(targetName, 'w');
    fprintf(fid, '%d\n', nPts);
    for j = 1:nSeg
%        fprintf(fid, '%f\n', datVals{j});
        fprintf(fid, '%0.16f\n', dt(j));
    end
    fprintf(fid, '%0.16f\n', dt(1));
    fclose(fid);
end