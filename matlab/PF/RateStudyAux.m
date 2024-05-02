function [valsRaw, matProcesseds, regScalars, regYs, PFs, storages, L1phi_delsig, maxLog, names_summary, names_summaryLatex, l_vals, lbl, studyName, studyMode] = ...
    RateStudyAux(studyMode, l_valsIn_s, isHyper, model_s, CZM_normalization4AT1_2, la_default_s, tauModel, includeHRAsymptotic, add_L1phi_delsig, plotData, optionNo)
max_x_size_4_plot = 6000000;
DMaxHyper = 0.99;

labsz = 25; % x, y label font size
lfs = 13;
printfig = 1;

% studyMode
% 1 loading rate
% 2 b - wave speed kept fixed - is set tauModel = lcoh2cd
% 3 b - taud kept fixed  OR the input
% 4 k (cd/cu ratio)

% tauModel = 
%           auto            -> if CZM_normalization4AT = 1 ->   tauModel = 'lcoh2cd'
%                              else                             tauModel = 'b2cd'

%           b2cd            tauScale = b/cd, cd = Phase field speed
%           lcoh2cd         tauScale = lcoh2cd/cd, cd = Phase field speed



if (nargin < 1)
    studyMode = 1;
end

% log10 of values in
% l_valsIn
if (nargin < 2)
    l_valsIn_s = {'-3', '-2.5', '-2', '-1.5', '-1', '-0.5', '0'};
    l_valsIn_s = {'-3', '-2', '-1', '0'};
end

if (nargin < 3)
    % the only valid options are 0 (PPF) and 1 (HPF)
    isHyper = 0;
end

if (nargin < 4)
    model_s = 'AT1';   % AT2, CZM-W, CZM-L
%    model_s = 'CZM-W';   % AT2, CZM-W, CZM-L
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
    add_L1phi_delsig = 0; % add energy of del sigma
end

if (nargin < 10)
    plotData = 1;
end

if (nargin < 11)
    optionNo = 0;
end

if (isnumeric(isHyper) == 0)
    if (strcmp(isHyper, 'E') == 1)
        isHyper = -1;
    elseif (strcmp(isHyper, 'P') == 1)
        isHyper = 0;
    elseif (strcmp(isHyper, 'H') == 1)
        isHyper = 1;
    end
end
if (isnumeric(studyMode) == 0)
    studyMode = str2num(studyMode);
end

if (studyMode == 2)
    tauModel = 'lcoh2cd';
end

CZM_modelName = 'Linear';

IMax = 2;
if (includeHRAsymptotic == 1)
    IMax = 3;
end

num = length(l_valsIn_s);
l_vals = zeros(num, 1);
vals = zeros(num, 1);
for i = 1:num
    l_vals(i) = str2num(l_valsIn_s{i});
    vals(i) = power(10.0, l_vals(i));
end
valsInv = 1.0 ./ vals;


isApproximate = 0;
la_s = '0';
l_cD2c_s = '0';
bTimesPiCZM_s = 'none';
df_s = '1';
vals_default0 = {isHyper, isApproximate, model_s, la_s, l_cD2c_s, df_s, CZM_normalization4AT1_2, bTimesPiCZM_s, CZM_modelName, tauModel};

bPrimeDefault = 1.0;
if (strcmp(model_s, 'CZM-W') == 1)
    bPrimeDefault = 2.0 / pi;
elseif (strcmp(model_s, 'CZM-L') == 1)
    bPrimeDefault = 1.0; % check
else
    if (CZM_normalization4AT1_2)
        if (strcmp(model_s, 'AT1') == 1)
            bPrimeDefault = 3.0 / 8.0;
        elseif (strcmp(model_s, 'AT2') == 1)
            bPrimeDefault = 27.0 / 256.0;
        end
    end
end

la_ss = cell(num, 1);
l_cD2c_ss = cell(num, 1);
bPrimes = cell(num, 1);
for i = 1:num
    la_ss{i} = la_default_s;
    l_cD2c_ss{i} = l_cD2c_s;
    bPrimes{i} = bPrimeDefault;
end


IStaticComputed = 1;
calculate_dvals_dInput = 0;
if (studyMode == 1) % loadingRate
    studyName = 'la';
    for i = 1:num
        la_ss{i} = l_valsIn_s{i};
    end
    lbl = '$$ \mathrm{log}({\bar{\dot{\epsilon}}}) $$';
    calculate_dvals_dInput = (strcmp(model_s, 'CZM-L') == 0);
elseif ((studyMode == 2) || (studyMode == 3))
    IStaticComputed = 0;
    if (studyMode == 2)
        studyName = 'lb_wsFixed';
    else
        studyName = 'lb_taudFixed';
    end
    for i = 1:num
        bPrimes{i} = vals(i) * bPrimeDefault;
    end
    lbl = '$$ \mathrm{\log}(b) $$';
elseif (studyMode == 4)
    studyName = 'lwsRatio';
    for i = 1:num
        l_cD2c_ss{i} = l_valsIn_s{i};
    end
    lbl = '$$ \mathrm{log}(k) $$';
end

if (IStaticComputed == 1)
    vals_default = vals_default0;
    vals_default{4} = 'none';
    vals_default{5} = '0';
    vals_default{1} = -1;
    bPrime = bPrimes{1};
    PFtmp = PF;
    [PF_IStatic, vals_out] = PFtmp.Initialize_Stage1(vals_default);
    PF_IStatic.bTimesPiCZM = -1;
    PF_IStatic.bPrime = bPrime;
    PF_IStatic = PF_IStatic.Compute();
end

nameBase = ['Study_option_', num2str(optionNo), '_mode_', num2str(studyMode), '_', studyName, '_H_', num2str(isHyper), '_model_', model_s, '_normCZM_', num2str(CZM_normalization4AT1_2)];
[status,message,messageid] = mkdir(nameBase);

% quasi-static solution

% I index: 0    exact solution, 1 quasi-static solution, 2 dynamic solution
% asymptotic


for i = 1:num
    vals_default1 = vals_default0;
    vals_default1{4} = la_ss{i};
    vals_default1{5} = l_cD2c_ss{i};
    bPrime = bPrimes{i};
    for I = 1:IMax
        vals_default = vals_default1;
        if (I == 2) % static
            vals_default{1} = -1;
            vals_default{4} = 'none';
        elseif (I == 3) % approximate high rate
            vals_default{2} = 1;
        end
        if ((I == 2) && IStaticComputed)
            PFs{i, I} = PF_IStatic;
        else
            PFtmp = PF;
            [PFtmp, vals_out] = PFtmp.Initialize_Stage1(vals_default);
            PFtmp.bTimesPiCZM = -1;
            PFtmp.bPrime = bPrime;
            PFtmp.DMaxHyper = DMaxHyper;
            PFtmp = PFtmp.Compute();
            PFs{i, I} = PFtmp;
        end
    end
end

colorNum = 1;
isModel = 0;
includeBlack = 0;
lc = getColors(colorNum, isModel, includeBlack);
num_vec = length(PFs{1, 1}.vecs_names);
vecNames = PFs{1, 1}.vecs_names;
num_scalars = length(PFs{1, 1}.scalar_names);
scalarNames = PFs{1, 1}.scalar_names;
scalar_names_latex = PFs{1, 1}.scalar_names_latex;
vecs_names_latex = PFs{1, 1}.vecs_names_latex;
xlab = vecs_names_latex{num_vec};

leg_names{1} = lbl;
for i = 1:length(l_valsIn_s);
    leg_names{i + 1} = l_valsIn_s{i};
end
I = 1;

if (plotData)
    close('all');
    for vi = 1:num_vec - 1
        vec_name = vecNames{vi};
        fnbase = [nameBase, '/', nameBase, '_A_eps_vs_', vec_name];
        fg = figure(vi);
        set(fg,'defaultLegendAutoUpdate', 'off');
        ylab = vecs_names_latex{vi};
        plot([nan], [nan], 'Color', 'w', 'LineStyle','none');
        hold on;

        xSize = 0;
        for ri = 1:num
            PFtmp = PFs{ri, I};
            x = PFtmp.epsilon_p;
            xSize = max(xSize, length(x));
            y = PFtmp.y_vecs{vi};
            plot(x, y, 'Color', lc{ri}, 'LineStyle', '-', 'LineWidth', 2);
            hold on;
        end
        lg = legend(leg_names, 'FontSize', lfs, 'Interpreter', 'latex');
        legend('boxoff');
 
        xh = get(gca, 'XLabel');
        set(xh, 'String', xlab, 'FontSize', labsz, 'VerticalAlignment','Top', 'Interpreter', 'latex');
        yh = get(gca, 'YLabel');
        set(yh, 'String', ylab, 'FontSize', labsz, 'VerticalAlignment','Bottom', 'Interpreter', 'latex');
    
        print('-dpng', [fnbase, '.png']);
        if (printfig)
            if (xSize < max_x_size_4_plot)
                savefig([fnbase, '.fig']);
            end
        end
    end
    close('all');
end


%%%%%%%%%%%%%%%%% plotting difference to zero loading rate ONLY for zero
%%%%%%%%%%%%%%%%% loading rate option

storages = cell(0);
L1phi_delsig = [];
if (calculate_dvals_dInput == 1)
    L1phi_delsig = zeros(num);
    xi = PF_IStatic.xi;
    offset = 2;
    if (xi == 1)
        offset = 3;
    end
    bBar = PF_IStatic.bPrime;
    for i = 1:num
        eps_vec = PFs{i, 1}.y_vecs{num_vec};
        sz = length(eps_vec);
        storage = zeros(num_vec, sz);
        for j = 1:sz
            storage(num_vec, j) = eps_vec(j);
        end
%        for j = 1:offset - 1
%            for vi = 1:num_vec - 1
%                PFtmp.y_vecs{vi}(j) = 0.0;
%            end
%        end
        for j = offset:sz
            eps = eps_vec(j);
            [vecStatic, fullDamage] = Compute_Elliptic_values_4_1_epsilon(xi, eps, bBar);
            for vi = 1:num_vec - 1
                y = PFs{i, 1}.y_vecs{vi}(j);
                y0 = vecStatic(vi);
                storage(vi, j) = y - y0;
            end
        end
        delsig = abs(storage(1, :));
        L1phi_delsig(i) = trapz(eps_vec, delsig);
        storages{i} = storage;
%        for vi = 1:num_vec - 1
%            tmpV = storage(vi, :);
%            PFtmp.y_vecs{vi} = tmpV;
%        end
%        PFDiffs{i} = PFtmp;
    end

    if (plotData)
        close('all');
        for vi = 1:num_vec - 1
            vec_name = vecNames{vi};
            fnbase = [nameBase, '/', nameBase, '_C_eps_vs_Der_', vec_name];
            fg = figure(vi + 20);
            set(fg,'defaultLegendAutoUpdate', 'off');
            ylab = vecs_names_latex{vi};
            plot([nan], [nan], 'Color', 'w', 'LineStyle','none');
            hold on;
    
            for ri = 1:num
                storage = storages{ri};
                x = storage(num_vec, :);
                y = storage(vi, :) * valsInv(ri);
%                PFtmp = PFDiffs{ri};
%                x = PFtmp.epsilon_p;
%                y = PFtmp.y_vecs{vi} * valsInv(ri);
                plot(x, y, 'Color', lc{ri}, 'LineStyle', '-', 'LineWidth', 2);
                hold on;
            end
            lg = legend(leg_names, 'FontSize', lfs, 'Interpreter', 'latex');
            legend('boxoff');
     
            xh = get(gca, 'XLabel');
            set(xh, 'String', xlab, 'FontSize', labsz, 'VerticalAlignment','Top', 'Interpreter', 'latex');
            yh = get(gca, 'YLabel');
            set(yh, 'String', ylab, 'FontSize', labsz, 'VerticalAlignment','Bottom', 'Interpreter', 'latex');
        
            print('-dpng', [fnbase, '.png']);
            if (printfig)
                savefig([fnbase, '.fig']);
            end
        end
        close('all');

    
    
        % now plotting diffs
        for vi = 1:num_vec - 1
            vec_name = vecNames{vi};
            fnbase = [nameBase, '/', nameBase, '_D_eps_vs_Diff_', vec_name];
            fg = figure(vi + 30);
            set(fg,'defaultLegendAutoUpdate', 'off');
            ylab = vecs_names_latex{vi};
            plot([nan], [nan], 'Color', 'w', 'LineStyle','none');
            hold on;
    
            xSize = 0;
            for ri = 1:num
%                PFtmp = PFDiffs{ri};
%                x = PFtmp.epsilon_p;
%                y = PFtmp.y_vecs{vi};
                storage = storages{ri};
                x = storage(num_vec, :);
                xSize = max(xSize, length(x));
                y = storage(vi, :);
                plot(x, y, 'Color', lc{ri}, 'LineStyle', '-', 'LineWidth', 2);
                hold on;
            end
            lg = legend(leg_names, 'FontSize', lfs, 'Interpreter', 'latex');
            legend('boxoff');
     
            xh = get(gca, 'XLabel');
            set(xh, 'String', xlab, 'FontSize', labsz, 'VerticalAlignment','Top', 'Interpreter', 'latex');
            yh = get(gca, 'YLabel');
            set(yh, 'String', ylab, 'FontSize', labsz, 'VerticalAlignment','Bottom', 'Interpreter', 'latex');
        
            print('-dpng', [fnbase, '.png']);
            if (printfig)
                if (xSize < max_x_size_4_plot)
                    savefig([fnbase, '.fig']);
                end
            end
        end
        close('all');
    end
end

%%%%% scalars and and convergence rate
has_phi_delsig = ((length(L1phi_delsig) > 0) && (add_L1phi_delsig));
nscalar = num_scalars + has_phi_delsig;

num_scalars_out = 2 * nscalar - 2;
if (has_phi_delsig)
    scalarNamesBK = scalarNames;
    scalar_names_latexBK = scalar_names_latex;
    scalarNames{1} = 'phi_del_sigma';
    scalar_names_latex{1} = '\phi(\Delta(\sigma))';
    for k = 1:num_scalars
        scalarNames{k + 1} = scalarNamesBK{k};
        scalar_names_latex{k + 1} = scalar_names_latexBK{k};
    end
end

for i = 1:nscalar - 2
    col(i) = i;
    isLog(i) = 1;
    names_summary{i} = ['log_', scalarNames{i}];
    names_summaryLatex{i} = ['\mathrm{log}(', scalar_names_latex{i}, ')'];
end

maxLog = i;
leg_reg_fields{1} = '$$ \mathrm{field} $$';
for i = 1:maxLog
    leg_reg_fields{i + 1} = ['$$ ', names_summaryLatex{i},' $$'];
end


st = i;
for j = 1:nscalar
    i = j + st;
    col(i) = j;
    isLog(i) = 0;
    names_summary{i} = [scalarNames{j}];
    names_summaryLatex{i} = [scalar_names_latex{j}];
end

valsRaw = cell(IMax, 1);
for I = 1:IMax
    for i = 1:num
        for j = 1:num_scalars
            mat(i, j + has_phi_delsig) = PFs{i, I}.scalars(j);
        end
    end
    if (has_phi_delsig)
        mat(i, 1) = 0;
        if (I == 2)
            for i = 1:num
                mat(i, 1) = L1phi_delsig(i);
            end
        end
    end
    valsRaw{I} = mat;
end

valsProcessed = cell(IMax, 1);
fnbase = [nameBase, '/', nameBase, '_B_Regressions_'];

for I = 1:IMax
    valsBase = valsRaw{1};
    valsCompared2 = valsRaw{I};
    if (I == 1)
        vals2Use = valsBase;
    else
        vals2Use = valsBase - valsCompared2;
    end
    matProcessed = zeros(num, num_scalars_out);
    for i = 1:num
        for j = 1:num_scalars_out
            coln = col(j);
            bLog = isLog(j);
            val = vals2Use(i, coln);
            if (bLog)
                vl = abs(val);
                if (vl > 1e-40)
                    val = log10(vl);
                else
                    val = nan;
                end
            end
            matProcessed(i, j) = val;
        end
    end
    matProcesseds{I} = matProcessed;
    fido = fopen([fnbase, '3_data_I_', num2str(I), '.csv'], 'w');
    fprintf(fido, 'names,%s', studyName);
    for j = 1:num_scalars_out
        fprintf(fido, ',%s', names_summary{j});
    end
    fprintf(fido, '\n');
    fprintf(fido, 'runNo,%s', lbl);
    for j = 1:num_scalars_out
        fprintf(fido, ',%s', names_summaryLatex{j});
    end
    for i = 1:num
        fprintf(fido, '\n%d,%g', i, l_vals(i));
        for j = 1:num_scalars_out
            fprintf(fido, ',%g', matProcessed(i, j));
        end
    end
    fclose(fido);
end

x = l_vals;
fid{1} = fopen([fnbase, '0_slope.csv'], 'w');
fid{2} = fopen([fnbase, '1_intercept.csv'], 'w');
fid{3} = fopen([fnbase, '2_r2.csv'], 'w');

for fi = 1:3
    fd = fid{fi};
    fprintf(fd, 'I');
    for j = 1:maxLog
        fprintf(fd, ',%s', scalarNames{j});
    end
end

for I = 1:IMax
    matProcessed = matProcesseds{I};
    for fi = 1:3
        fd = fid{fi};
        fprintf(fd, '\n%d', I);
    end
    for j = 1:maxLog
        y = matProcessed(:, j);
        % Perform linear regression
        p = polyfit(x, y, 1);
        slope = p(1);
        intercept = p(2);
        % Calculate R-squared value
        yfit = polyval(p, x);
        yresid = y - yfit;
        SSresid = sum(yresid.^2);
        SStotal = (length(y)-1) * var(y);
        rsq = 1 - SSresid/SStotal;
        regs = [slope, intercept, rsq];

        for fi = 1:3
            fd = fid{fi};
            fprintf(fd, ',%g', regs(fi));
        end
        regScalars{I}{j} = regs;
        regYs{I}{j} = yfit;
    end
end
for fi = 1:3
    fclose(fid{fi});
end

if (plotData)
    close('all');
    for I = 1:IMax
        matProcessed = matProcesseds{I};
        fnbase = [nameBase, '/', nameBase, '_B_regression_', num2str(I)];
        fg = figure(10 + I);
        set(fg,'defaultLegendAutoUpdate', 'off');
        plot([nan], [nan], 'Color', 'w', 'LineStyle','none');
        hold on;

        for j = 1:maxLog
            y = matProcessed(:, j);
            plot(x, y, 'Color', lc{j}, 'LineStyle', '-', 'LineWidth', 2);
            hold on;
        end
        lg = legend(leg_reg_fields, 'FontSize', lfs, 'Interpreter', 'latex');
        legend('boxoff');

        for j = 1:maxLog
            y = regYs{I}{j};
            plot(x, y, 'Color', lc{j}, 'LineStyle', '--', 'LineWidth', 2);
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
    end
    close('all');
end

