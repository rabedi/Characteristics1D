function [l10epsDots, l10sMs, l10phids] = get_PF0D_Responses_EPH_PF_loading_Rates(model_s, bVal, bMode, la_del, la_min, la_max, ...
G, E, sigmac)
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
    bVal = 1.0;
end
if (nargin < 3)
    bMode = 2;
end
if (nargin < 4)
    la_del = 0.25;
end
if (nargin < 5)
    la_min = -2;
end
if (nargin < 6)
    la_max = 6;
end
if (nargin < 7)
    G = 1.0;
end
if (nargin < 8)
    E = 1.0;
end
if (nargin < 9)
    sigmac = 1.0;
end
l10epsDots = la_min:la_del:la_max;
sz = length(l10epsDots);
for i = 1:sz
    la_ss{i} = num2str(l10epsDots(i));
end

isHyper = -1;
la_s = '0';
[sMEPF, phidEPF, PFtmp_EPF] = get_PF0D_Response(model_s, isHyper, bVal, bMode, la_s, G, E, sigmac);
l10sMEPF = log10(sMEPF);
l10phidEPF = log10(phidEPF);
l10sMs = l10sMEPF * ones(3, sz);
l10phids = l10phidEPF * ones(3, sz);

for isHyper = 0:1
    row = isHyper + 2;
    fprintf(1, '\nisHyper: %d', isHyper);
    for i = 1:sz
        la_s = la_ss{i};
        fprintf(1, ' %s', la_s);
        [sMEPF, phidEPF, PFtmp] = get_PF0D_Response(model_s, isHyper, bVal, bMode, la_s, G, E, sigmac);
        l10sMs(row, i) = log10(sMEPF);
        l10phids(row, i) = log10(phidEPF);
    end
end
