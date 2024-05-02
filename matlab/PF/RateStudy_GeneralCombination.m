function RateStudy_GeneralCombination(l_valsIn_s, studyModes, isHypers, model_ss, la_default_ss, ...
    optionNo, CZM_normalization4AT1_2, tauModel, ...
    legendEntriesLatex, header, linesClrStyleMarkers, regressionChangeMode, ...
    includeHRAsymptotic, add_L1phi_delsig, plotData, multiply_breaker_ab)
labsz = 25; % x, y label font size
lfs = 13;
printfig = 1;

% log10 of values in
% l_valsIn
if (nargin < 1)
    l_valsIn_s = {'-3', '-2.5', '-2', '-1.5', '-1', '-0.5', '0'};
    l_valsIn_s = {'-3', '-2', '-1'};
%    l_valsIn_s = {'-4', '-3', '-2', '-1'};
%    l_valsIn_s = {'-4', '-3', '-2'};
end

if (nargin < 2)
    studyModes = {'1'};
end

if (nargin < 3)
    isHypers = {'P'};
end

if (nargin < 4)
    model_ss = {'AT1', 'AT2', 'CZM-W'}; % 'CZM-L
    model_ss = {'AT1', 'CZM-W'}; % 'CZM-L
end

if (nargin < 5)
    la_default_ss = {'0'};
end


if (nargin < 6)
    optionNo = 0;
end

if (nargin < 7)
    CZM_normalization4AT1_2 = 0;
end

if (nargin < 8)
    tauModel = 'auto';
end

if (nargin < 9)
    legendEntriesLatex = cell(0);
end

if (nargin < 10)
    header = 'auto';
end
if (nargin < 11)
    linesClrStyleMarkers = [];
end

if (nargin < 12)
    regressionChangeMode = -1;
    % -1 -> decide based on number of varying breakers
    % 0 -> line style is changed to --
    % 1 -> marker is off
end

if (nargin < 13)
    includeHRAsymptotic = 0;
end

if (nargin < 14)
    add_L1phi_delsig = 1;
end

if (nargin < 15)
    plotData = 1;
end

if (nargin < 16)
    % this is used to multiply field breaker a and b; see its use below 
    multiply_breaker_ab = [];
end
%%%%%%%%%%%%%%%
IMax = 2;
if (includeHRAsymptotic == 1)
    IMax = 3;
end

% Initialization
isModel = 0;
isHPE = 0;

num = length(l_valsIn_s);
l_vals = zeros(num, 1);
vals = zeros(num, 1);
for i = 1:num
    l_vals(i) = str2num(l_valsIn_s{i});
    vals(i) = power(10.0, l_vals(i));
end

breakerVals{1} = studyModes;
breakerVals{2} = isHypers;
breakerVals{3} = model_ss;
breakerVals{4} = la_default_ss;
breakerValsBK = breakerVals;

szm = length(multiply_breaker_ab);
for bi = 1:4
    szs(bi) = length(breakerVals{bi});
end
if (szm > 1)
    a = multiply_breaker_ab(1);
    b = multiply_breaker_ab(2);
    veca = breakerVals{a};
    vecb = breakerVals{b};
    cntr = 0;
    for ia = 1:szs(a)
        for ib = 1:szs(b)
            cntr = cntr + 1;
            breakerVals{a}{cntr} = veca{ia};
            breakerVals{b}{cntr} = vecb{ib};
        end
    end
    szs(a) = cntr;
    szs(b) = cntr;
end 

breakerNames = {'study', 'H', 'model', 'la'};

if (szs(1) == 0)
    breakerVals{1}{1} = '1';
    szs(1) = 1;
end
if (szs(2) == 0)
    breakerVals{2}{1} = '0';
    szs(2) = 1;
end

if (szs(3) == 0)
    breakerVals{3}{1} = 'AT1';
    szs(3) = 1;
end

if (szs(4) == 0)
    breakerVals{4}{1} = '0';
    szs(4) = 1;
end

a = find(szs > 1);
num_Varying = length(a);
runName = '';
for bi = 1:4
    if (szs(bi) > 1)
        continue;
    end
    if (~isempty(runName))
        runName = [runName, '_'];
    end
    nm = breakerVals{bi}{1};
    runName = [runName, breakerNames{bi}, '_', nm];
end

numCombos = max(szs);
for bi = 1:4
    sz = szs(bi);
    st = sz + 1;
    for ci = st:numCombos
        ciiBack = mod(ci, sz) + 1;
        breakerVals{bi}{ci} = breakerVals{bi}{ciiBack};
    end
end

b_legendEntriesLatex = (length(legendEntriesLatex) > 0);
b_header = (strcmp(header, 'auto') == 0);
b_linesClrStyleMarkers = (length(linesClrStyleMarkers) > 0);
b_regressionChangeMode = (regressionChangeMode < 0);

if (num_Varying == 1)
    breaker = a(1);
    isModel = (breaker == 3);
    isHPE = (breaker == 2);

    if (b_header == 0)
        header = breakerNames{breaker};
    end
    if (b_legendEntriesLatex == 0)
        for ci = 1:numCombos
            legendEntriesLatex{ci} = breakerVals{breaker}{ci};
        end
    end
    if ((b_linesClrStyleMarkers == 0) || (b_regressionChangeMode == 0))
        regressionChangeMode = 0; % line style is changed to --
        for ci = 1:numCombos
            linesClrStyleMarkers{ci} = [ci, 1, 1];
        end
    end
elseif (num_Varying == 2)
    breaker0 = a(1);
    breaker1 = a(2);

    if (b_header == 0)
        header = [breakerNames{breaker0}, ',', breakerNames{breaker1}];
    end
    if (b_legendEntriesLatex == 0)
        for ci = 1:numCombos
            legendEntriesLatex{ci} = [breakerVals{breaker0}{ci}, ',', breakerVals{breaker1}{ci}];
        end
    end
    if ((b_linesClrStyleMarkers == 0) || (b_regressionChangeMode == 0))
        regressionChangeMode = 1; % marker is removed
        for ci = 1:numCombos
            vl0 = breakerVals{breaker0}{ci};
            vl1 = breakerVals{breaker1}{ci};
            cit = 1;
            for cit = 1:length(breakerValsBK{breaker0})
                if (strcmp(breakerValsBK{breaker0}{cit}, vl0) == 1)
                    break;
                end
            end
            cit0 = cit;
            cit = ci;
            for cit = 1:length(breakerValsBK{breaker1})
                if (strcmp(breakerValsBK{breaker1}{cit}, vl1) == 1)
                    break;
                end
            end
            cit1 = cit;
            
            linesClrStyleMarkers{ci} = [cit0, 1, cit1 + 1];
        end
    end
else
    if (b_header == 0)
        header = 'header';
    end
    if (b_legendEntriesLatex == 0)
        fprintf(1, 'legendEntriesLatex not provided\n');
        pause;
    end
    if ((b_linesClrStyleMarkers == 0) || (b_regressionChangeMode == 0))
        fprintf(1, 'b_linesClrStyleMarkers or b_regressionChangeMode not provided\n');
        pause;
    end
end

%%%%%%%%%%%%%%% line colors, .... 
includeBlack = (isModel || isHPE);
lc = getColors(1, isModel, includeBlack);

lstyles = {'-', '--', '-.', ':'};
lths = [2, 2, 2, 3];
lMarkers = {'none', 'o', 's', 'd', 'p', '<', '>'};


for ci = 1:numCombos
    studyMode = breakerVals{1}{ci};
    isHyper = breakerVals{2}{ci};
    model_s = breakerVals{3}{ci};
    la_default_s = breakerVals{4}{ci};
    [valsRaw{ci}, matProcesseds{ci}, regScalars{ci}, regYs{ci}, PFs{ci}, storages{ci}, L1phi_delsig{ci}, maxLog, names_summary, names_summaryLatex, l_vals, lbl, studyName, studyMode] = ...
    RateStudyAux(studyMode, l_valsIn_s, isHyper, model_s, CZM_normalization4AT1_2, la_default_s, tauModel, includeHRAsymptotic, add_L1phi_delsig, plotData, optionNo);
end

nameBase = ['Study_option_', num2str(optionNo), '_', runName, '_normCZM_', num2str(CZM_normalization4AT1_2)];
[status,message,messageid] = mkdir(nameBase);


%%% printing and plotting things
fnbase = [nameBase, '/', nameBase, '_B_Regressions_'];
num_scalars_out = length(names_summary);

% part 1, scalar values for all runs
for I = 1:IMax
    fido = fopen([fnbase, '3_data_I_', num2str(I), '.csv'], 'w');
    fprintf(fido, 'runNo,breakerNo,breakerSymbol,val');
    for j = 1:num_scalars_out
        fprintf(fido, ',%s', names_summary{j});
    end
    fprintf(fido, '\n');
    fprintf(fido, 'runNo,breakerNo,breakerSymbol,%s', lbl);
    for j = 1:num_scalars_out
        fprintf(fido, ',%s', names_summaryLatex{j});
    end
    for i = 1:num
        for ci = 1:numCombos
            matProcessed = matProcesseds{ci}{I};
            symbol = legendEntriesLatex{ci};
            symbol = strrep(symbol, ',', '-');        
            fprintf(fido, '\n%d,%d,%s,%g', i, ci, symbol, l_vals(i));
            for j = 1:num_scalars_out
                fprintf(fido, ',%g', matProcessed(i, j));
            end
        end
    end
    fclose(fido);
end
% regression info
fnbase = [nameBase, '/', nameBase, '_B_Regressions_'];
fid{1} = fopen([fnbase, '0_slope.csv'], 'w');
fid{2} = fopen([fnbase, '1_intercept.csv'], 'w');
fid{3} = fopen([fnbase, '2_r2.csv'], 'w');

scalarNames = names_summary;
for fi = 1:3
    fd = fid{fi};
    fprintf(fd, 'I,breakerNo,breakerSymbol');
    for j = 1:maxLog
        fprintf(fd, ',%s', scalarNames{j});
    end
end

IMax = 2;
if (includeHRAsymptotic == 1)
    IMax = 3;
end

x = l_vals;
for I = 1:IMax
    for ci = 1:numCombos
        symbol = legendEntriesLatex{ci};
        symbol = strrep(symbol, ',', '-');        
        for fi = 1:3
            fd = fid{fi};
            fprintf(fd, '\n%d,%d,%s', I, ci, symbol);
        end
        regss = regScalars{ci}{I};
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

legendEntry{1} = header;
for ci = 1:numCombos
    legendEntry{ci + 1} = legendEntriesLatex{ci};
end

numFields = length(names_summary);
if (plotData)
    close('all');
    for I = 1:IMax
        for j = 1:numFields
            isLog = (j <= maxLog);
            fg = figure(j);
            fnbase = [nameBase, '/', nameBase, '_B_regressionOrData_I_', num2str(I), '_fld_', num2str(j), '_', names_summary{j}];
            set(fg,'defaultLegendAutoUpdate', 'off');
            plot([nan], [nan], 'Color', 'w', 'LineStyle','none');
            hold on;
            for ci = 1:numCombos
                matProcessed = matProcesseds{ci}{I};
                lProp = linesClrStyleMarkers{ci};
                lc_n = lProp(1);
                lc_style_n = lProp(2);
                lc_marker_n = lProp(3);
                lcv = lc{lc_n};
                lstylev = lstyles{lc_style_n};
                lth = lths(lc_style_n);
                mrkr = lMarkers{lc_marker_n};

                y = matProcessed(:, j);
                plot(x, y, 'Color', lcv, 'LineStyle', lstylev, 'LineWidth', lth, 'Marker', mrkr, 'MarkerEdgeColor', lcv);
                hold on;
            end
            lg = legend(legendEntry, 'FontSize', lfs, 'Interpreter', 'latex');
            legend('boxoff');

            if (isLog)
                for ci = 1:numCombos
                    y = regYs{ci}{I}{j};
    
                    lProp = linesClrStyleMarkers{ci};
                    lc_n = lProp(1);
                    lc_style_n = lProp(2);
                    lc_marker_n = lProp(3);
                    lcv = lc{lc_n};
                    
                    if (regressionChangeMode == 0)
                        lc_style_n = 2;
                    elseif (regressionChangeMode == 1)
                        lc_marker_n = 1;
                    end
                    lstylev = lstyles{lc_style_n};
                    lth = lths(lc_style_n);
                    mrkr = lMarkers{lc_marker_n};
    
                    plot(x, y, 'Color', lcv, 'LineStyle', lstylev, 'LineWidth', lth, 'Marker', mrkr, 'MarkerEdgeColor', lcv);
                    hold on;
                end
            end
 
            xh = get(gca, 'XLabel');
            set(xh, 'String', lbl, 'FontSize', labsz, 'VerticalAlignment','Top', 'Interpreter', 'latex');
            yh = get(gca, 'YLabel');
            ylab = ['$$ ', names_summaryLatex{j}, ' $$'];
            set(yh, 'String', ylab, 'FontSize', labsz, 'VerticalAlignment','Bottom', 'Interpreter', 'latex');
    
            print('-dpng', [fnbase, '.png']);
            if (printfig)
                savefig([fnbase, '.fig']);
            end
            close('all');
        end
    end
end

