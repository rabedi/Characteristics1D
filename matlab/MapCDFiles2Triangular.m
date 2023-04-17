function MapCDFiles2Triangular(alpha, beta, log2NumSegments, numRealizations)
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
    numRealizations = 500;
end
nSeg = 2^log2NumSegments;
nSegSer = num2str(nSeg);
nPts = nSeg + 1;
nPtsSer = num2str(nPts);

rootFolder = ['../InhomogeneousFiles/CD/al', alpha, 'b', beta];
targetFolder = [rootFolder, '_np', nPtsSer];
%targetFolder = [rootFolder, '_np', nSegSer];

fido = fopen(['summaryCD_al', alpha, '_b', beta, '.csv'], 'w');
dd2 = 0.5;

for i = 1:numRealizations
    ii = i - 1;
    targetName = [targetFolder, '/initial_values_', num2str(ii), '.txt'];
    fid = fopen(targetName, 'r');
    if (fid < 0)
        return;
    end
    nPts = fscanf(fid, '%d', 1);
    dat = zeros(nPts, 1);
    for j = 1:nSeg
        dat(j) = fscanf(fid, '%f', 1);
    end
    for j = 1:nSeg
        v = dat(j);
        triVal = SN2Tri(v, dd2);
        fprintf(fido, '%f', triVal);
        if (j < nSeg)
            fprintf(fido, ',');
        else
            fprintf(fido, '\n');
        end
    end
    fclose(fid);
end
fclose(fido);