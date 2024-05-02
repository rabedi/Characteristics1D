l_valsIn_s = {'-3', '-2', '-1'};
% studyMode
% 1 loading rate
% 2 b - wave speed kept fixed - is set tauModel = lcoh2cd
% 3 b - taud kept fixed  OR the input
% 4 k (cd/cu ratio)

studyModes = {'1'};
isHypers = {'P'};

model_ss = {'AT1', 'AT2', 'CZM-W'}; % 'CZM-L
model_ss = {'AT1', 'CZM-W'}; % 'CZM-L

la_default_ss = {'0'};

optionNo = 0;

CZM_normalization4AT1_2 = 0;

tauModel = 'auto';

legendEntriesLatex = cell(0);
header = 'auto';
linesClrStyleMarkers = [];
regressionChangeMode = -1;

includeHRAsymptotic = 0;
add_L1phi_delsig = 1;
plotData = 1;
multiply_breaker_ab = []; % multiply modelss by is Hyper

if (0) % plot P,H & model at the same time
    optionNo = 10;
    CZM_normalization4AT1_2 = 1;
    studyModes = {'2'}; % b changing

    model_ss = {'AT1', 'AT2', 'CZM-W'}; % 'CZM-L
    model_ss = {'AT1', 'CZM-W'}; % 'CZM-L
    isHypers = {'P', 'H'};
    % study -> 1, isHypers -> 2, model -> 3, la_defaults -> 4
    multiply_breaker_ab = [2, 3]; % multiply modelss by is Hyper
    RateStudy_GeneralCombination(l_valsIn_s, studyModes, isHypers, model_ss, la_default_ss, ...
        optionNo, CZM_normalization4AT1_2, tauModel, ...
        legendEntriesLatex, header, linesClrStyleMarkers, regressionChangeMode, ...
        includeHRAsymptotic, add_L1phi_delsig, plotData, multiply_breaker_ab);
end


if (0) % plot P,H & model at the same time
    optionNo = 20;
    CZM_normalization4AT1_2 = 1;
    studyModes = {'2'}; % b changing

    model_ss = {'AT1', 'AT2', 'CZM-W'}; % 'CZM-L
 %   model_ss = {'AT1', 'CZM-W'}; % 'CZM-L
%    isHypers = {'E', 'P', 'H'};
    isHypers = {'E'};
    RateStudy_GeneralCombination(l_valsIn_s, studyModes, isHypers, model_ss, la_default_ss, ...
        optionNo, CZM_normalization4AT1_2, tauModel, ...
        legendEntriesLatex, header, linesClrStyleMarkers, regressionChangeMode, ...
        includeHRAsymptotic, add_L1phi_delsig, plotData, multiply_breaker_ab);
end

if (0) % plot P,H & model at the same time
%    model_ssBase = {'AT1', 'AT2', 'CZM-W'}; % 'CZM-L
    model_ssBase = {'AT1'}; % 'CZM-L
    l_valsIn_s = {'-3', '-2.5', '-2', '-1.5', '-1', '-0.5', '0'};
    la_default_ss = {'0'}; % 0 too
%    l_valsIn_s = {'-3', '-2.5', '-2', '-1.5', '-1', '-0.5', '0'};
    la_default_ss = {'2'};
    l_valsIn_s = {'-4', '-3.5', '-3', '-2.5', '-2', '-1.5', '-1', '-0.5', '0'};
%    la_default_ss = {'3'};
%    l_valsIn_s = {'-5', '-4.5', '-4', '-3.5', '-3', '-2.5', '-2', '-1.5', '-1', '-0.5', '0'};
%    la_default_ss = {'4'};
    for j = 1:length(model_ssBase)
        model_ss = cell(1);
        model_s = model_ssBase{j};
        model_ss{1} = model_s;
        if (strcmp(model_s, 'AT1') == 1)
            offset = 1;
        elseif (strcmp(model_s, 'AT2') == 1)
            offset = 2;
        elseif (strcmp(model_s, 'CZM-W') == 1)
            offset = 3;
        elseif (strcmp(model_s, 'CZM-L') == 1)
            offset = 4;
        end
        optionNo = 30 + offset;
        CZM_normalization4AT1_2 = 0;
        studyModes = {'2'}; % b changing
        isHypers = {'E', 'P', 'H'};
        RateStudy_GeneralCombination(l_valsIn_s, studyModes, isHypers, model_ss, la_default_ss, ...
            optionNo, CZM_normalization4AT1_2, tauModel, ...
            legendEntriesLatex, header, linesClrStyleMarkers, regressionChangeMode, ...
            includeHRAsymptotic, add_L1phi_delsig, plotData, multiply_breaker_ab);
    end
end



if (1) % plot P,H & model at the same time
    optionNo = 100;
    CZM_normalization4AT1_2 = 1;
    studyModes = {'1'}; % a changing

    model_ss = {'AT1', 'AT2', 'CZM-W'}; % 'CZM-L
 %   model_ss = {'AT1', 'CZM-W'}; % 'CZM-L
%    isHypers = {'E', 'P', 'H'};
    isHypers = {'P', 'H'};
    multiply_breaker_ab = [2, 3]; % multiply modelss by is Hyper
    RateStudy_GeneralCombination(l_valsIn_s, studyModes, isHypers, model_ss, la_default_ss, ...
        optionNo, CZM_normalization4AT1_2, tauModel, ...
        legendEntriesLatex, header, linesClrStyleMarkers, regressionChangeMode, ...
        includeHRAsymptotic, add_L1phi_delsig, plotData, multiply_breaker_ab);
end
