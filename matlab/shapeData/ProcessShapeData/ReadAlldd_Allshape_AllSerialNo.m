function ReadAlldd_Allshape_AllSerialNo(serNoMax)
if nargin < 1
    serNoMax = 4999;
end
ddStrs = {'0.1', '0.3', '0.5', '0.7', '0.9'};
n_dd = length(ddStrs);
for di = 1:n_dd
    Read1dd_Allshape_AllSerialNo(ddStrs{di}, serNoMax)
end
