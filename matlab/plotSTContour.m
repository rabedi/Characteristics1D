function plotSTContour(X, T, vals, numSpatialSegmentsFinal, hasFracture)
[nt, nx] = size(X);
if (hasFracture)
    numSpatialSegments = nx / 2;
else
    numSpatialSegments = nx - 1;
end
step = floor(numSpatialSegments / numSpatialSegmentsFinal);
if (step < 2)
    contourf(X, T, vals, 'EdgeColor', 'none');
    return;
end
xinds = [];
if (hasFracture)
    cntr = 1;
    xinds(cntr) = 1;
    for j = step:step:numSpatialSegments - 1
        cntr = cntr + 1;
        pos = j * 2;
        xinds(cntr) = pos;
        cntr = cntr + 1;
        xinds(cntr) = pos + 1;
    end
    cntr = cntr + 1;
    xinds(cntr) = nx;
else
    cntr = 1;
    xinds(cntr) = 1;
    for j = step:step:numSpatialSegments - 1
        cntr = cntr + 1;
        pos = 1 + j;
        xinds(cntr) = pos;
    end
    cntr = cntr + 1;
    xinds(cntr) = nx;
end    
XRed = X(:,xinds);
TRed = T(:, xinds);
valsRed = vals(:, xinds);
contourf(XRed, TRed, valsRed, 'EdgeColor', 'none');
