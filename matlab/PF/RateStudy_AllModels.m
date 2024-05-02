function RateStudy_AllModels(studyMode, l_valsIn_s, isHyper, model_ss, CZM_normalization4AT1_2, la_default_s, tauModel, includeHRAsymptotic, add_L1phi_delsig, plotData)
labsz = 25; % x, y label font size
lfs = 13;
printfig = 1;

if (nargin < 1)
    studyMode = 2;
end

% log10 of values in
% l_valsIn
if (nargin < 2)
    l_valsIn_s = {'-3', '-2.5', '-2', '-1.5', '-1', '-0.5', '0'};
    l_valsIn_s = {'-3', '-2', '-1'};
%    l_valsIn_s = {'-4', '-3', '-2', '-1'};
    l_valsIn_s = {'-4', '-3', '-2'};
end

if (nargin < 3)
    % the only valid options are 0 (PPF) and 1 (HPF)
    isHyper = 0;
end

if (nargin < 4)
    model_ss = {'AT1', 'AT2', 'CZM-W'}; % 'CZM-L
end

if (nargin < 5)
    CZM_normalization4AT1_2 = 0;
end

if (nargin < 6)
    la_default_s = '0';
end

if (nargin < 7)
    tauModel = 'auto';
end

if (nargin < 8)
    includeHRAsymptotic = 0;
end
if (nargin < 9)
    add_L1phi_delsig = 1; % add energy of del sigma
end

if (nargin < 10)
    plotData = 1;
end

num_models = length(model_ss);
for mi = 1:num_models
    model_s = model_ss{mi};
    [valsRaw{mi}, matProcesseds{mi}, regScalars{mi}, regYs{mi}, PFs{mi}, storages{mi}, L1phi_delsig{mi}, maxLog, names_summary, names_summaryLatex, l_vals, lbl, studyName, studyMode] = ...
    RateStudyAux(studyMode, l_valsIn_s, isHyper, model_s, CZM_normalization4AT1_2, la_default_s, tauModel, includeHRAsymptotic, add_L1phi_delsig, plotData);
end

nameBase = ['Study', num2str(studyMode), '_', studyName, '_H_', num2str(isHyper), '_normCZM_', num2str(CZM_normalization4AT1_2)];
[status,message,messageid] = mkdir(nameBase);

fnbase = [nameBase, '/', nameBase, '_B_Regressions_'];
fid{1} = fopen([fnbase, '0_slope.csv'], 'w');
fid{2} = fopen([fnbase, '1_intercept.csv'], 'w');
fid{3} = fopen([fnbase, '2_r2.csv'], 'w');

scalarNames = names_summary;
for fi = 1:3
    fd = fid{fi};
    fprintf(fd, 'I,model_s');
    for j = 1:maxLog
        fprintf(fd, ',%s', scalarNames{j});
    end
end


x = l_vals;
for I = 1:IMax
    for mi = 1:num_models
        model_s = model_ss{mi};
        for fi = 1:3
            fd = fid{fi};
            fprintf(fd, '\n%d,%s', I, model_s);
        end
        regss = regScalars{mi}{I};
        for j = 1:maxLog
            regs = regss{j};
            for fi = 1:3
                fd = fid{fi};
                fprintf(fd, ',%g', regs(fi));
            end
        end
    end
end
for fi = 1:3
    fclose(fid{fi});
end

isModel = 1;
includeBlack = 1;
lc = getColors(1, isModel, includeBlack);

legendEntry{1} = 'model';
for mi = 1:num_models
    model_s = model_ss{mi};
    legendEntry{mi + 1} = model_s;
end

if (plotData)
    close('all');
    for I = 1:IMax
        for j = 1:maxLog
            fg = figure(j);
            fnbase = [nameBase, '/', nameBase, '_B_regression_', num2str(I), '_fld_', num2str(j)];
            set(fg,'defaultLegendAutoUpdate', 'off');
            plot([nan], [nan], 'Color', 'w', 'LineStyle','none');
            hold on;
            for mi = 1:num_models
                model_s = model_ss{mi};
                matProcessed = matProcesseds{mi}{I};

                y = matProcessed(:, j);
                plot(x, y, 'Color', lc{mi}, 'LineStyle', '-', 'LineWidth', 2);
                hold on;
            end
            lg = legend(legendEntry, 'FontSize', lfs, 'Interpreter', 'latex');
            legend('boxoff');

            for mi = 1:num_models
                y = regYs{mi}{I}{j};
                plot(x, y, 'Color', lc{mi}, 'LineStyle', '--', 'LineWidth', 2);
                hold on;
            end
 
            xh = get(gca, 'XLabel');
            set(xh, 'String', lbl, 'FontSize', labsz, 'VerticalAlignment','Top', 'Interpreter', 'latex');
            yh = get(gca, 'YLabel');
            set(yh, 'String', 'ylabel', 'FontSize', labsz, 'VerticalAlignment','Bottom', 'Interpreter', 'latex');
    
            print('-dpng', [fnbase, '.png']);
            if (printfig)
                savefig([fnbase, '.fig']);
            end
            close('all');
        end
    end
end

