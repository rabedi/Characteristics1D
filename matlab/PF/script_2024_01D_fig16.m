model_s = 'CZMW';
la_del = 0.25;
la_min = -2;
la_max = 6.5;

isHyper = -1;
bVal = 1.0;
%       0: b2l = 1 (AT1, AT2, l = lc)
%       1: actual b provided
%       2: b2lc Relative to sM = 1 b provided (e.g. bVal = 1 for AT1 ->
%       b2lc = 3/8) ...
%       3: b2lc: actual value is provided
bMode = 2.0;
la_s = '0';
G = 1;
E = 1;
sigmac = 1;



%%%%% actual paper values
%bVal = 200e-9;
%bMode = 1;

bVal = 1.0/32.0;
%bVal = 1.0/4.0;
bMode = 2;

G = 83.13;
E = 610e9;
sigmac = 1e9;

[sM, phid, PFtmp] = get_PF0D_Response(model_s, isHyper, bVal, bMode, la_s, G, E, sigmac);
sM
phid

[l10epsDots, l10sMs, l10phids] = get_PF0D_Responses_EPH_PF_loading_Rates(model_s, bVal, bMode, la_del, la_min, la_max, ...
G, E, sigmac);
save('fig_16_l10epsDots.txt', 'l10epsDots', '-ascii');
save('fig_16_l10sMs.txt', 'l10sMs', '-ascii');
save('fig_16_l10phids.txt', 'l10phids', '-ascii');

figure(1);
plot(l10epsDots, l10sMs);
print('-dpng', 'fig16_l10epsDots_l10sMs.png');
figure(2);
plot(l10epsDots, l10phids);
print('-dpng', 'fig16_l10epsDots_l10phids.png');


