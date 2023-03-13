function plot_contours(displacementScale, serialNum, versionNum, subdomainNos, hasFracture, plotDeformed, numSpatialSegmentsFinal, root, plotOpenedUpRing_uv, equalAxis)
fs = 18;
if nargin < 1
    displacementScale = 0.2;
%    displacementScale = 1.0;
end
if nargin < 2
    serialNum = 0;
end
if nargin < 3
    versionNum = -1;
    versionNum = 14;
end
if nargin < 4
    subdomainNos = [0];
end
if nargin < 5
    hasFracture = -1;
end
if nargin < 6
    plotFractures = 1;
end
if nargin < 7
    plotDeformed = 1;
end
if nargin < 8
    numSpatialSegmentsFinal = -1; % inactive -> no resolution reduction
%    numSpatialSegmentsFinal = 512;
%    numSpatialSegmentsFinal = 256;
%    numSpatialSegmentsFinal = 32;
end
if nargin < 9
    root = '../';
end
if nargin < 10
    % if 1 plots v' = vtheta + ax and correspondign v (that is the field we are directly solving for)
    %    0 plots vtheta which is already printed in files
    % this applies to corresponding displacement too
    plotOpenedUpRing_uv = 1;
    plotOpenedUpRing_uv = 0;
end
if nargin < 11
    equalAxis = 1;
end
if nargin < 12
    plotFig = 0;
end
folderName = [root, 'run'];
if (versionNum >= 0)
    folderName = [folderName, '_V_', num2str(versionNum)];
end
folderName = [folderName, '_', num2str(serialNum), '/'];
fn = [folderName, '_visualization_timesteps.txt'];
fid = fopen(fn, 'r');
buf = fscanf(fid, '%s', 1);
nt = fscanf(fid, '%d', 1);
buf = fscanf(fid, '%s', 1);
timStep = fscanf(fid, '%f', 1);
for ni = 1:nt
    ts(ni) = fscanf(fid, '%f', 1);
end
fclose(fid);

figureBaseName = [folderName, '_plot_sds_'];
pF = 0;
for sii = 1:length(subdomainNos)
    si = subdomainNos(sii);
    figureBaseName = [figureBaseName, num2str(sii)];
    sd = ['_sd_', num2str(si), '_'];
    fbase = [folderName, sd, 'visualization_'];
    s = load([fbase, 's.txt'], '-ascii');
    v = load([fbase, 'v.txt'], '-ascii');
    u = load([fbase, 'u.txt'], '-ascii');
    [numT, numX] = size(s);
    fn = [fbase, 'paras.txt'];
    fid = fopen(fn, 'r');
    buf = fscanf(fid, '%s', 1);
    a = fscanf(fid, '%f', 1);
    buf = fscanf(fid, '%s', 1);
    nx = fscanf(fid, '%d', 1);
    if (nx ~= numX)
        fprintf(1, 'nx (%d) != numX (%d)\n', nx, numX);
        pause
    end
    if (nt ~= numT)
%        fprintf(1, 'nt (%d) != numT (%d)\n', nt, numT);
%        pause
        ts = [];
        for j = 1:numT
            ts(j) = (j - 1) * timStep;
        end
    end
    bf = fscanf(fid, '%s', 1);
    for j = 1:numX
        xs(j) = fscanf(fid, '%f', 1);
    end
    L = max(xs) - min(xs);

    fclose(fid);

    hF = hasFracture;
    if (hasFracture == -1)
        hF = ((xs(3) - xs(2)) < 1e-7 * L);
    end
    pF = hF &&  plotFractures;
    fn = [fbase, 'D.txt'];
    D = [];
    DdelUMax = [];
    DdelU = [];
    if (pF)
        fd = fopen(fn, 'r');
        if (fd < 0)
            pF = 0;
        else
            D = load([fbase, 'D.txt'], '-ascii');
            DdelUMax = load([fbase, 'DdelUMax.txt'], '-ascii');
            DdelU = load([fbase, 'DdelU.txt'], '-ascii');
            fclose(fd);
        end
    end
    
    [X, T] = meshgrid(xs, ts);
    if (plotOpenedUpRing_uv)
        aX = a * X;
        aXT = T .* aX;
        u = u + aXT;
        v = v + aX;
    end
    cntr = 1;
    figure(cntr);
    plotSTContour(X, T, s, numSpatialSegmentsFinal, hF);
    hold on;
    cntr = cntr + 1;
    figure(cntr);
    plotSTContour(X, T, v, numSpatialSegmentsFinal, hF);
    hold on;
    cntr = cntr + 1;
    figure(cntr);
    plotSTContour(X, T, u, numSpatialSegmentsFinal, hF);
    hold on;
    cntr = cntr + 1;
    if (pF)
        figure(cntr);
        plotSTContour(X, T, D, numSpatialSegmentsFinal, hF);
        hold on;
        cntr = cntr + 1;

        figure(cntr);
        plotSTContour(X, T, DdelUMax, numSpatialSegmentsFinal, hF);
        hold on;
        cntr = cntr + 1;

        figure(cntr);
        plotSTContour(X, T, DdelU, numSpatialSegmentsFinal, hF);
        hold on;
        cntr = cntr + 1;
    end
    
    if (plotDeformed)
        Xdeformed = X + displacementScale * u;
        figure(cntr);
        plotSTContour(Xdeformed, T, s, numSpatialSegmentsFinal, hF);
        hold on;
        cntr = cntr + 1;
        figure(cntr);
        plotSTContour(Xdeformed, T, v, numSpatialSegmentsFinal, hF);
        hold on;
        cntr = cntr + 1;
    end
end
figureBaseName = [figureBaseName, '_'];

fldNames = {'s', 'v', 'u'};
fldTitleNames = fldNames;
cntr = length(fldNames) + 1;
if (pF)
    fldNames{cntr} = 'D';
    fldTitleNames{cntr} = 'D';
    cntr = cntr + 1;
    fldNames{cntr} = 'DdelUMax';
    fldTitleNames{cntr} = 'DdelUMax';
    cntr = cntr + 1;
    fldNames{cntr} = 'DdelU';
    fldTitleNames{cntr} = 'DdelU';
    cntr = cntr + 1;
end
if (plotDeformed)
    fldNames{cntr} = 's';
    fldTitleNames{cntr} = 'sDeformed';
    cntr = cntr + 1;
    fldNames{cntr} = 'v';
    fldTitleNames{cntr} = 'vDeformed';
    cntr = cntr + 1;
end
sz = length(fldTitleNames);

for fi = 1:sz
    fn = [figureBaseName, fldTitleNames{fi}];
    figure(fi);
    colorbar;
    hx = xlabel('$$ x $$', 'interpreter', 'latex', 'FontSize', fs);
    set(hx, 'VerticalAlignment','Top');
    hy = ylabel('$$ t $$', 'interpreter', 'latex', 'FontSize', fs);
    set(hy, 'VerticalAlignment','Bottom');
    if (equalAxis)
        axis('equal');
    end
    print('-dpng', [fn, '.png']);
    if (plotFig)
        savefig([fn, '.fig']);
    end
end
fclose('all');