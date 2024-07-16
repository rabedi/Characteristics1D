fclose('all');
close('all');
clear all;
ldf = 0;
add0D = 1;

labsz = 25; % x, y label font size
lfs = 11;

%l10saG = [-1, 0, 1, 1.5, 2, 2.5, 3, 3.5];
%l10saG = [0, 1];
%l10saG = [-1, -0.5];
l10saG = [6];
l10saG = [-1, -0.5];
l10saG = [0.5, 1];
%l10saG = [0, 0.5];
l10sa = l10saG - 3;
epsDot = '$$ \mathrm{log}_{10}({\dot{\bar{\epsilon}}}_{{\circ}H}) ';
nla = length(l10sa);

cntr = 0;
hyperSym = {'EPF', 'PPF', 'HPF'};
clrs = {...
    [0 1 0], ... %green
    [0 135/255 0], ... %green
    [0, 0.5, 1], ... 
    [0	0	1], ... % blue 'b'  ;%    [0	0	0.5], ... % dark blue    
    [1	102/255	0], ... %orange'
    [1	0	0] ...
 };
%    [0.5	0.25	0], ... % brown
%    [0.5	0.5	0.5], ...%'dark_gray2'
%    [0 135/255 0], ... %green


szPDE = 3;
for isHyper = -1:1
    hs = [hyperSym{isHyper + 2}, ', '];
    for i = 1:nla
        la = l10sa(i);
        cntr = cntr + 1;
        [eps_vec{cntr}, sig_vec{cntr}, epsM(cntr), sigM(cntr), epsF(cntr), phiF(cntr)] = ReadGiangData_NonDimIO(la, ldf, isHyper);
        leg{cntr} = [hs, epsDot, '=', num2str(la), ' $$'];
        lav(cntr) = la;
    end
end
cntr = cntr + 1;
leg{cntr} = '0D,\ PPF\ \&\ HPF';
cntr = cntr + 1;
leg{cntr} = '0D,\ EPF';

if (add0D)
    isApproximate = 0;
    model_s = 'AT1';
    la_s = '0';
    l_cD2c_s = '0';
    df_s = '0';
    CZM_normalization4AT1_2 = 0;
    bTimesPiCZM_s = 'none';
    CZM_modelName = 'linear';
    tauModel = 'auto';
    vals_default0 = {isHyper, isApproximate, model_s, la_s, l_cD2c_s, df_s, CZM_normalization4AT1_2, bTimesPiCZM_s, CZM_modelName, tauModel};

    cntr = 0;
    for isHyper = -1:1
        for i = 1:nla
            cntr = cntr + 1;
            la = l10sa(i);
            vals_default = vals_default0;
            vals_default{1} = isHyper;
            vals_default{4} = num2str(la);
            PFtmp = PF;
            [PFtmp, vals_out] = PFtmp.Initialize_Stage1(vals_default);
            PFtmp = PFtmp.Compute();
            PFs{cntr} = PFtmp;
        end
    end
end


%colorNum = 1;
%isModel = 0;
%includeBlack = 0;
%lc = getColors(colorNum, isModel, includeBlack);

close('all');
fnbase = 'compare0D1D_PH';
fg = figure(1);
set(fg,'defaultLegendAutoUpdate', 'off');
%plot([nan], [nan], 'Color', 'w', 'LineStyle','none');
%hold on;


sz = nla * szPDE;
for i = 1:sz
    la = lav(i);
    x = eps_vec{i};
    y = sig_vec{i};
    plot(x, y, 'Color', clrs{i}, 'LineStyle', '-', 'LineWidth', 2);
    hold on;
end
if (add0D)
    plot([nan], [nan], 'Color', [0.5	0.5	0.5], 'LineStyle','--', 'LineWidth', 3);
    hold on;
    plot([nan], [nan], 'Color', 'k', 'LineStyle','-.', 'LineWidth', 3);
    hold on;
end

lg = legend(leg, 'FontSize', lfs, 'Interpreter', 'latex');
legend('boxoff');

if (add0D)
    for i = 1:sz
        ls = '--';
        lc = clrs{i};
        PFtmp = PFs{i};
        if ((i > 1) && (PFtmp.isHyper == -1))
            continue;
        elseif (i == 1)
            ls = '-.';
            lc = 'k';
        end
        x = PFtmp.epsilon_p;
        y = PFtmp.sigma_p;
        plot(x, y, 'Color', lc, 'LineStyle', ls, 'LineWidth', 3);
        hold on;
    end
end
xlim([0, 2.5]);
%xlim([0, 1.75]);
ylim([-0.07, 0.67]);
xlab = '$$ \bar{\epsilon} $$';
ylab = '$$ \bar{\sigma} $$';
xh = get(gca, 'XLabel');
set(xh, 'String', xlab, 'FontSize', labsz, 'VerticalAlignment','Top', 'Interpreter', 'latex');
yh = get(gca, 'YLabel');
set(yh, 'String', ylab, 'FontSize', labsz, 'VerticalAlignment','Bottom', 'Interpreter', 'latex');

print('-dpng', [fnbase, '.png']);
savefig([fnbase, '.fig']);
%close('all');
