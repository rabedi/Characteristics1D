function [sM, phid, PFtmp] = get_PF0D_Response(model_s, isHyper, bVal, bMode, la_s, G, E, sigmac)
% bMode: 
%       0: b2l = 1 (AT1, AT2, l = lc)
%       1: actual b provided
%       2: b2lc Relative to sM = 1 b provided (e.g. bVal = 1 for AT1 ->
%       b2lc = 3/8) ...
%       3: b2lc: actual value is provided
if (nargin < 1)
    model_s = 'AT1';
end
if (nargin < 2)
    isHyper = 1;
end
if (nargin < 3)
    bVal = 1.0;
end
if (nargin < 4)
    bMode = 2;
end
if (nargin < 5)
    la_s = '0';
end
if (nargin < 6)
    G = 1.0;
end
if (nargin < 7)
    E = 1.0;
end
if (nargin < 8)
    sigmac = 1.0;
end
isAT1 = strcmp(model_s, 'AT1');
isAT2 = strcmp(model_s, 'AT2');
isCZMW = strcmp(model_s, 'CZMW');
isCZML = strcmp(model_s, 'CZML');
b2lcBase = 1.0;
if (isAT1)
    b2lcBase = 3.0/8.0;
elseif (isAT2)
    b2lcBase = 27.0/256.0;
elseif (isCZMW)
    b2lcBase = 1.0 / pi;
elseif (isCZML)
    b2lcBase = 0.25;
end

isCZM = isCZML || isCZMW;
b2lc = 1.0;
b = 1.0;
lc = G * E / sigmac / sigmac;
if (bMode == 0) % modeling AT1, AT2, for lTilde = l
    b = bVal;
    b2lc = 1.0; % it's reall b2l
elseif (bMode == 1) % actual val given
    b = bVal;
    b2lc = bVal / lc;
elseif (bMode == 3) % b2lc ratio is given
    b2lc = bVal;
    b = b2lc * lc;
elseif (bMode == 2) % like to above but b2lc is given relative to default vals
    b2lc = b2lcBase * bVal;
    b = b2lc * lc;
end
isApproximate = 0;
l_cD2c_s = 'none';
df_s = 'none';
CZM_normalization4AT = 0;
bTimesPiCZM_s = '1.0';
CZM_modelName = 'Linear';
%tauModel = 'auto';
%df = power(10, -0.5);
%df_s = num2str(df);
tauModel = 'b2cd';
tauModel = 'lcoh2cd';

vals_default = {isHyper, isApproximate, model_s, la_s, l_cD2c_s, df_s, CZM_normalization4AT, bTimesPiCZM_s, CZM_modelName, tauModel};
PFtmp = PF;
[PFtmp, vals_out] = PFtmp.Initialize_Stage1(vals_default);
PFtmp.bPrime = b2lc;
PFtmp.bTimesPiCZM = -1.0;
PFtmp = PFtmp.Compute();
sM = PFtmp.sigmap_Max;
phid = PFtmp.phi;





