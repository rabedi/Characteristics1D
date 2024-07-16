fclose('all');
close('all');
useMarker = 0;
plotProblemAtLow_df_HPF = 1; % generates Hlod_dfProb_all ... show cases the problem with HPF, low df
%option = 30;
lwdd = 2;
%isHyper = -1;
isHyper = 1;
isApproximate = 0;
tauModelIn = 'b2cd';
% AT1, AT2, CZM-W, CZM-L are the options
model_ss = {'AT1', 'CZM-W'};
model_ss = {'AT1'};
model_ss = {'CZM-W'};
num_model_s = length(model_ss);
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



PFtmp = cell(0);

close('all');
for mi = 1:num_model_s
    model_s = model_ss{mi};
    vals_default = {isHyper, isApproximate, model_s, la_s, l_cD2c_s, df_s, CZM_normalization4AT, bTimesPiCZM_s, CZM_modelName, tauModelIn};
    PFtmp{mi} = PF;
    [PFtmp{mi}, vals_out] = PFtmp{mi}.Initialize_Stage1(vals_default);
    PFtmp{mi} = PFtmp{mi}.Compute();
%    figure(1);
%    plot(PFtmp.epsilon_p, PFtmp.sigma_p);
%    figure(2);
%    plot(PFtmp.epsilon_p, PFtmp.D_vec);
    
    xH{mi} = PFtmp{mi}.epsilon_p;
    yH{mi}{1} = PFtmp{mi}.sigma_p;
    yH{mi}{2} = PFtmp{mi}.D_vec;
end

PFtmp2 = cell(0);
if plotProblemAtLow_df_HPF
for mi = 1:num_model_s
    model_s = model_ss{mi};
    isHyper = -1;
    vals_default = {isHyper, isApproximate, model_s, la_s, l_cD2c_s, df_s, CZM_normalization4AT, bTimesPiCZM_s, CZM_modelName, tauModelIn};
    
    PFtmp2{mi} = PF;
    [PFtmp2{mi}, vals_out] = PFtmp2{mi}.Initialize_Stage1(vals_default);
    PFtmp2{mi} = PFtmp2{mi}.Compute();
%    figure(11);
%    plot(PFtmp2.epsilon_p, PFtmp2.sigma_p);
%    figure(12);
%    plot(PFtmp2.epsilon_p, PFtmp2.D_vec);
    
    xE{mi} = PFtmp2{mi}.epsilon_p;
    yE{mi}{1} = PFtmp2{mi}.sigma_p;
    yE{mi}{2} = PFtmp2{mi}.D_vec;
end

yLabel{1} = '$$ {\bar{\sigma}}_M $$';
yLabel{2} = '$$ d $$';


lfs = 14;
labsz = 25; % x, y label font size
ch = [0.5	0.5	0.5];
cp = [0 0 0];
xLabel = '$$ {\bar{\dot{\epsilon}}}_0 $$';
leg_names = {'HPF', 'EPF'};
fnbase = 'Hlod_dfProb';
for mi = 1:num_model_s
    fn = [fnbase, '_', model_ss{mi}];
    for i = 1:2
        is = num2str(i);
        figure(3 + i);
        plot(xH{mi}, yH{mi}{i}, '-', 'Color', ch, 'LineWidth', 2);
        hold on;
        plot(xE{mi}, yE{mi}{i}, '-', 'Color', cp, 'LineWidth', 2);
        hold on;
    
        lg = legend(leg_names, 'FontSize', lfs);
        legend('boxoff');
        
        xh = get(gca, 'XLabel');
        set(xh, 'String', xLabel, 'FontSize', labsz, 'VerticalAlignment','Top', 'Interpreter', 'latex');
        yh = get(gca, 'YLabel');
        set(yh, 'String', yLabel{i}, 'FontSize', labsz, 'VerticalAlignment','Bottom', 'Interpreter', 'latex');
        xlim([0, 6]);
        
        print('-dpng', [fn, is, '.png']);
        savefig([fn, is, '.fig']);
    end
end




%%%
if (num_model_s == 1)
    fn = [fnbase, '_', model_ss{1}];
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
    
        plot(xH{1}, yH{1}{axi}, '-', 'Color', clr, 'LineWidth', 2);
        hold on;
        plot(xE{1}, yE{1}{axi}, '--', 'Color', clr, 'LineWidth', 2);
        hold on;
        ylim([0, 1]);
        yh = get(gca, 'YLabel');
        set(yh, 'String', yLabel{axi}, 'FontSize', labsz, 'Color', clr, 'VerticalAlignment', vert, 'Interpreter', 'latex');
    end    
        
    %lg = legend(leg_names, 'FontSize', lfs);
    %legend('boxoff');
        
    xh = get(gca, 'XLabel');
    set(xh, 'String', xLabel, 'FontSize', labsz, 'VerticalAlignment','Top', 'Interpreter', 'latex');
    if (strcmp(model_ss{1}, 'CZM-W') == 1)
        xlim([0, 2.2]);
    else
        xlim([0, 6]);
    end
    print('-dpng', [fn, '_all.png']);
    savefig([fn, '_all.fig']);

else
    leg_namesBK = {'HPF', 'EPF'};
    n_leg_names = 2;
    leg_names = cell(0);
    mrkrs = {'s', 'o', 'd', 'p'};
    if ((num_model_s == 2) && (useMarker == 0))
        lstyles = {'-', '-.', '--', ':'};
        leg_names = {[model_ss{1}, ', HPF'], [model_ss{1}, ', EPF'], [model_ss{2}, ', HPF'], [model_ss{2}, ', EPF']};
        l_markers = {'none', 'none', 'none', 'none'};
    else
        cntr = 0;
        for mi = 1:num_model_s
            model_s = model_ss{mi};
            for j = 1:n_leg_names
                cntr = cntr + 1;
                leg_names{cntr} = [model_s, ' ,', leg_namesBK{j}];
                l_markers{cntr} = mrkrs{mi};
                lst = '-';
                if (j == 2)
                    lst = '--';
                end
                lstyles{cntr} = lst;
            end
        end
    end
    n_lines = length(lstyles);

    fn = fnbase;
    close('all');
    fg = figure(8);
    set(fg,'defaultLegendAutoUpdate', 'off');
    for li = 1:n_lines
        ls = lstyles{li};
        lw = 2;
        if (ls == ':')
            lw = lwdd;
        end
        plot([nan], [nan], 'Color', [0.5 0.5 0.5], 'LineStyle', ls, 'LineWidth', lw, 'Marker', l_markers{li});
        hold on;
    end
    lg = legend(leg_names, 'FontSize', lfs);
    legend('boxoff');
    
    ax = gca;
    for axi = 1:2
        cntr = 0;
        for mi = 1:num_model_s
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
        
            cntr = cntr + 1;
            ls = lstyles{cntr};
            lw = 2;
            if (ls == ':')
                lw = lwdd;
            end
            plot(xH{mi}, yH{mi}{axi}, 'LineStyle', ls, 'Color', clr, 'LineWidth', lw, 'Marker', l_markers{cntr});
            hold on;
            cntr = cntr + 1;
            ls = lstyles{cntr};
            lw = 2;
            if (ls == ':')
                lw = lwdd;
            end
            plot(xE{mi}, yE{mi}{axi}, 'LineStyle', ls, 'Color', clr, 'LineWidth', lw, 'Marker', l_markers{cntr});
            hold on;
            ylim([0, 1]);
            yh = get(gca, 'YLabel');
            set(yh, 'String', yLabel{axi}, 'FontSize', labsz, 'Color', clr, 'VerticalAlignment', vert, 'Interpreter', 'latex');
        end
    end
        
    %lg = legend(leg_names, 'FontSize', lfs);
    %legend('boxoff');
        
    xh = get(gca, 'XLabel');
    set(xh, 'String', xLabel, 'FontSize', labsz, 'VerticalAlignment','Top', 'Interpreter', 'latex');
    if (strcmp(model_ss{1}, 'CZM-W') == 1)
        xlim([0, 2.2]);
    else
        xlim([0, 6]);
    end
print('-dpng', [fn, '_all.png']);
savefig([fn, '_all.fig']);


end


end