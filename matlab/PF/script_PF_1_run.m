plotProblemAtLow_df_HPF = 1; % generates Hlod_dfProb_all ... show cases the problem with HPF, low df
%option = 30;
%isHyper = -1;
isHyper = 0;
isApproximate = 0;
tauModelIn = 'b2cd';
% AT1, AT2, CZM-W, CZM-L are the options
model_s = 'AT1';
model_s = 'CZM-L';
la_s = 'none';
l_cD2c_s = 'none';
df_s = 'none';
CZM_normalization4AT = -1;
bTimesPiCZM_s = '1.0';
CZM_modelName = 'Linear';

%la_s = '-3.5';
% Damage becomes negative
la_s = '6';
df_s = '100';
isHyper = 1;

% eps_f < eps_f EPH
la_s = '-2';
df_s = '0.01';





close('all');
vals_default = {isHyper, isApproximate, model_s, la_s, l_cD2c_s, df_s, CZM_normalization4AT, bTimesPiCZM_s, CZM_modelName, tauModelIn};
PFtmp = PF;
[PFtmp, vals_out] = PFtmp.Initialize_Stage1(vals_default);
PFtmp = PFtmp.Compute();
figure(1);
plot(PFtmp.epsilon_p, PFtmp.sigma_p);
figure(2);
plot(PFtmp.epsilon_p, PFtmp.D_vec);

xH = PFtmp.epsilon_p;
yH{1} = PFtmp.sigma_p;
yH{2} = PFtmp.D_vec;



if plotProblemAtLow_df_HPF
isHyper = -1;
vals_default = {isHyper, isApproximate, model_s, la_s, l_cD2c_s, df_s, CZM_normalization4AT, bTimesPiCZM_s, CZM_modelName, tauModelIn};

PFtmp2 = PF;
[PFtmp2, vals_out] = PFtmp2.Initialize_Stage1(vals_default);
PFtmp2 = PFtmp2.Compute();
figure(11);
plot(PFtmp2.epsilon_p, PFtmp2.sigma_p);
figure(12);
plot(PFtmp2.epsilon_p, PFtmp2.D_vec);

xE = PFtmp2.epsilon_p;
yE{1} = PFtmp2.sigma_p;
yE{2} = PFtmp2.D_vec;

yLabel{1} = '$$ {\bar{\sigma}}_M $$';
yLabel{2} = '$$ d $$';


lfs = 14;
labsz = 25; % x, y label font size
ch = [0.5	0.5	0.5];
cp = [0 0 0];
xLabel = '$$ {\bar{\dot{\epsilon}}}_0 $$';
leg_names = {'HPF', 'EPF'};
fnbase = 'Hlod_dfProb';
for i = 1:2
    is = num2str(i);
    figure(3 + i);
    plot(xH, yH{i}, '-', 'Color', ch, 'LineWidth', 2);
    hold on;
    plot(xE, yE{i}, '-', 'Color', cp, 'LineWidth', 2);
    hold on;

    lg = legend(leg_names, 'FontSize', lfs);
    legend('boxoff');
    
    xh = get(gca, 'XLabel');
    set(xh, 'String', xLabel, 'FontSize', labsz, 'VerticalAlignment','Top', 'Interpreter', 'latex');
    yh = get(gca, 'YLabel');
    set(yh, 'String', yLabel{i}, 'FontSize', labsz, 'VerticalAlignment','Bottom', 'Interpreter', 'latex');
    xlim([0, 6]);
    
    print('-dpng', [fnbase, is, '.png']);
    savefig([fnbase, is, '.fig']);
end





%%%
close('all');
fg = figure(8);
set(fg,'defaultLegendAutoUpdate', 'off');
plot([nan], [nan], 'Color', [0.5 0.5 0.5], 'LineStyle', '-', 'LineWidth', 2);
hold on;
plot([nan], [nan], 'Color', [0.5 0.5 0.5], 'LineStyle', '--', 'LineWidth', 2);
hold on;
lg = legend(leg_names, 'FontSize', lfs);
legend('boxoff');

ax = gca;
for axi = 1:2
    clr = 'b';
    vert = 'Bottom';
    if (axi == 2)
        yyaxis right;
        clr = 'r';
        vert = 'Top';
    else
        yyaxis left;
    end
    ax.YColor = clr;

    plot(xH, yH{axi}, '-', 'Color', clr, 'LineWidth', 2);
    hold on;
    plot(xE, yE{axi}, '--', 'Color', clr, 'LineWidth', 2);
    hold on;
    ylim([0, 1]);
    yh = get(gca, 'YLabel');
    set(yh, 'String', yLabel{axi}, 'FontSize', labsz, 'Color', clr, 'VerticalAlignment', vert, 'Interpreter', 'latex');
end    
    
%lg = legend(leg_names, 'FontSize', lfs);
%legend('boxoff');
    
xh = get(gca, 'XLabel');
set(xh, 'String', xLabel, 'FontSize', labsz, 'VerticalAlignment','Top', 'Interpreter', 'latex');
if (strcmp(model_s, 'CZM-W') == 1)
    xlim([0, 2.2]);
else
    xlim([0, 6]);
end
print('-dpng', [fnbase, '_all.png']);
savefig([fnbase, '_all.fig']);
end