function [X, Y, xreg, yreg, R2, shape, scale, N] = loadWeibullPlotData4InhomogFiles(shapeStr, ddStr, llcStr, b_print, fld)
if nargin < 1
    shapeStr = '1.5';
end
if nargin < 2
    ddStr = '0.5';
end
if nargin < 3
    llcStr = '-1.5';
end
if nargin < 4
    b_print = 1;
end
if nargin < 5
    fld = 2;
end
fldStr = num2str(fld);
fileNWOExt = ['Weibull/fld', fldStr, '/shape', shapeStr, '_dd', ddStr, '_llc', llcStr, '_fld', fldStr, '_Weibull'];
[X, Y, xreg, yreg, R2, shape, scale, N] = loadWeibullPlotData(fileNWOExt);
if (b_print)
    plot(X, Y, '-k');
    hold on;
    plot(xreg, yreg, '-b');
    R2
    shape
    scale
    N
end