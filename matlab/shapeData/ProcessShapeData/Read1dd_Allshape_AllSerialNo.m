function Read1dd_Allshape_AllSerialNo(ddStr, serNoMax)
if nargin < 1
    ddStr = '0.5';
end
if nargin < 2
    serNoMax = 4999;
end
fprintf(1, 'dd = %s\n', ddStr);
shapeStrs = {'1', '1.5', '2', '3', '4'};
n_shape = length(shapeStrs);
for si = 1:n_shape
    [dat_llc_indexed, serNo_llc_indexed, cntrs_llc_indexed] = Read1dd_1shape_AllSerialNo(shapeStrs{si}, ddStr, serNoMax);
end
