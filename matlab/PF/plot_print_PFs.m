function plot_print_PFs(option, isHyper, isApproximate, model_s, la_s, l_cD2c_s, df_s, CZM_normalization4AT, bTimesPiCZM_s, CZM_modelName, tauModel, drawPlots)
labsz = 25; % x, y label font size
% does not add HRA, PPR, except EPF
additionalLines4_50 = 0;

%option 0 to 10 are for comparing various PF formulations
%option 0: AT1
%option 1: AT1, AT2
%option 2: CZM-Wu, CZM-Lorentz
%option 5: AT1, AT2, CZM-Wu, CZM-Lorentz
%option 6: AT1, AT2, CZM-Wu, CZM-Lorentz - no CZM normalization
%option 7: AT1, CZM-Wu, CZM-Lorentz


%option 10: CZM-Wu: for all TSR models (linear, ...)

%option 20: CZM-Wu: linear model, for different b' -> time scale = b/c_u ->
%not realistic
%option 21: time scale = lcoh / c_u

%%%%%%%% LOADING RATE EFFECT
% option 30-31,32: loading rate in legend (H/P) x (model type = AT1, ...) [cD2c
% and df may change too]
% 30 makes the plots: both approximate and exact
% 31 -> only exact
% 32 -> only approximate
%   IF CZM_normalization4AT == 0 -> 5 is added to option outputs [35, 36, 37]

%%%%%%%% DAMPING FACTOR EFFECT
% options 40-42:
% 40 for plottig, adds approximate lines if isApproximate ~= 0
% 41 -> only exact
% 42 -> only approximate
%   IF CZM_normalization4AT == 0 -> 5 is added to option outputs [45, 46, 47]
             
%%%%%%%% WAVE SPEED FACTOR
% options 50-52:
% 50 for plottig, adds approximate lines if isApproximate ~= 0
% 51 -> only exact
% 52 -> only approximate
%   IF CZM_normalization4AT == 0 -> 5 is added to option outputs [55, 56, 57]

%%%%%%%%%%%%%%%%%%%
% ap is the loading rate:
% if > 0, it's real, if < 0, for H and P PF, it computes the asymptotic
% high loading rate solution

% isHyper
%       -1: Elliptic
%        0: Parabolic
%        1: Hyperbolic

if nargin < 1
    option = 20; %5, 2, 10, 20
%    option = 21;
%    option = 5; % model comparison - CZM normalization
    option = 6; % model comparison - CZM normalization
    option = 30;
%    option = 40;
    option = 40;
    option = 50;
    option = 5;
    option = 6;
end


if nargin < 2
    isHyper = -1;
%    isHyper = 0;
%    isHyper = 1;
end

if nargin < 3
    % for dynamic solutions approximate = 0 -> plots the exact 0D solution,
    % if 1 it plots the High Loading rate limit
    isApproximate = 2;          % use to plot HRA plots
    % isApproximate = 1;        % HRA approximate part is added to initial
    % part (that is, it's not the HRA only solution)
end

if nargin < 4
    % AT1, AT2, CZM-W, CZM-L are the options
    model_s = 'AT1';
%    model_s = 'AT2';
%    model_s = 'CZM-W';
end


if nargin < 5
    la_s = 'none';
    la_s = '0';
    la_s = '2';
    la_s = '1';
%    la_s = '-2';
    la_s = 'none';
end

if nargin < 6
    l_cD2c_s = 'none';
end

if nargin < 7
    df_s = 'none';
end

if nargin < 8
    % 0 -> don't normalize
    % 1 -> ATs have the same sigma max as CZMs
    % -1 -> does not change the value in the class
    CZM_normalization4AT = -1;
    CZM_normalization4AT = 0;
    CZM_normalization4AT = 1;
end

if nargin < 9
    bTimesPiCZM_s = '1.0';
end

if nargin < 10
    CZM_modelName = 'Linear';
end

% tauModel = 
%           auto            -> if CZM_normalization4AT = 1 ->   tauModel = 'lcoh2cd'
%                              else                             tauModel = 'b2cd'

%           b2cd            tauScale = b/cd, cd = Phase field speed
%           lcoh2cd         tauScale = lcoh2cd/cd, cd = Phase field speed
if nargin < 11
    tauModel = 'auto';
    tauModel = 'b2cd';
%    tauModel = 'lcoh2cd';
end

if nargin < 12
    drawPlots = 1;
end
plot_D_xAxis = 0; % plots D vs. fields


vals_default = {isHyper, isApproximate, model_s, la_s, l_cD2c_s, df_s, CZM_normalization4AT, bTimesPiCZM_s, CZM_modelName, tauModel};
name_vals = {'isHyper', 'Approx', 'model', 'la', 'cD2c', 'df', 'normalize2CZM', 'bTimesPi', 'CZMmodelName', 'tauModel'};
nameHP = {'PPF', 'HPF'};
HR_name = 'HR';
for i = 1:2
    nameHP_HR{i} = ['HRA-', nameHP{i}];
end

vals_out = cell(0);
name_out = '';

% there is an option to specify these
lnums = cell(0);
lstys = cell(0);

lnum_extras = cell(0);
lstys_extras = cell(0);
names_extras = cell(0);

ap = -1;
la = -101;
if (strcmp(la_s, 'none') == 0)
    la = str2num(la_s);
    ap = power(10, la);
end


header = 'header';

% each curve is distniguished by a color
useColor4Lines = 1;
names = cell(0);
PFs = cell(0);
if (option < 10)
    header = 'model';
    CZM_normalization4AT1_2 = CZM_normalization4AT;
    if (CZM_normalization4AT1_2 == -1)
        CZM_normalization4AT1_2 = 1;
    end
    xis = [];
    omegaCZMs = [];
    if (option == 0) % AT1
        xis = [1];
        omegaCZMs = [0];
        names = {'AT1'};
        model_ss = {'AT1'};
        CZM_normalization4AT1_2 = 0;
    elseif (option == 1) T
        xis = [1, 0];
        omegaCZMs = [0, 0];
        names = {'AT1', 'AT2'};
        model_ss = {'AT1', 'AT2'};
        CZM_normalization4AT1_2 = 0;
    elseif (option == 2) % CZM-Wu, CZM-Lorentz
        xis = [2, 1];
        omegaCZMs = [1, 1];
        names = {'CZM-W', 'CZM-L'};
        model_ss = {'CZM-W', 'CZM-L'};
    elseif ((option == 5) || (option == 6)) % AT1, AT2, CZM-Wu, CZM-Lorentz
        xis = [1, 0, 2, 1];
        omegaCZMs = [0, 0, 1, 1];
        names = {'AT2', 'AT1', 'CZM-W', 'CZM-L'};
        model_ss = {'AT2', 'AT1', 'CZM-W', 'CZM-L'};
        if (option == 5)
            CZM_normalization4AT1_2 = 1;
        else
            CZM_normalization4AT1_2 = 0;
        end
    elseif (option == 7) % AT1, CZM-Wu, CZM-Lorentz
        xis = [1, 2, 1];
        omegaCZMs = [0, 1, 1];
        names = {'AT1', 'CZM-W', 'CZM-L'};
        model_ss = {'AT1', 'CZM-W', 'CZM-L'};
    end
    num = length(xis);
    PFs = cell(num, 1);
    for i = 1:num
        PFtmp = PF;

        vals_default = {isHyper, isApproximate, model_ss{i}, la_s, l_cD2c_s, df_s, CZM_normalization4AT1_2, bTimesPiCZM_s, CZM_modelName, tauModel};
        [PFtmp, vals_out{i}] = PFtmp.Initialize_Stage1(vals_default);
%         PFtmp.ap = ap;
%         PFtmp.isHyper = isHyper;
%         PFtmp.CZM_normalization4AT1_2 = CZM_normalization4AT1_2;
%         PFtmp.xi = xis(i);
%         PFtmp.omegaCZM = omegaCZMs(i);
%         PFtmp.CZM_model_name = 'Linear';
%         PFtmp.bPrime = -1;
         PFs{i} = PFtmp;
    end
elseif (option == 10)
    names = {'Linear', 'Bilinear', 'Exponential', 'Hyperbolic', 'Concrete'};
    num = length(names);
    for i = 1:num
        PFtmp = PF;
        PFtmp.ap = ap;
        PFtmp.isHyper = isHyper;
        PFtmp.CZM_normalization4AT1_2 = 1;
        PFtmp.xi = 2;
        PFtmp.omegaCZM = 1;
        PFtmp.CZM_model_name = names{i};
        PFtmp.bPrime = -1;
        PFs{i} = PFtmp;
    end
elseif ((option == 20) || (option == 21))

    tauModel = 'b2cd';
    if (option == 21)
        tauModel = 'lcoh2cd';
    end
    isCZM = (strncmp(model_s, 'CZM', 3) == 1);
    isAT1 = (strcmp(model_s, 'AT1') == 1);
    if (isHyper == -1)
        names = {'2', '1', '0.5'}; % bPrime * pi
    else
        if (isCZM)
            names = {'2', '1', '0.5', '0.25', '0.125', '0.0625', '0.03125'}; % bPrime * pi
            names = {'2', '1', '0.5', '0.25', '0.125', '0.0625', '0.03125', '0.015625'}; % bPrime * pi
        else
            names = {'1', '0.5', '0.25', '0.125', '0.0625', '0.03125', '0.015625'}; % bPrime * pi
        end
    end
%    names = {'0.0625', '0.03125'}; % bPrime * pi
    num = length(names);
    %model_s = 'CZM-W';
    CZM_modelName = 'Linear';
    for i = 1:num
        PFtmp = PF;
        PFtmp.ap = ap;
        PFtmp.isHyper = isHyper;
        PFtmp.CZM_normalization4AT1_2 = 1;
        PFtmp.xi = 2;
        PFtmp.omegaCZM = 1;
        PFtmp.CZM_model_name = 'Linear';
%        PFtmp.bPrime = str2num(names{i}) / pi;
%        PFtmp.bTimesPiCZM = -1;
        bTimesPiCZM_s = names{i};
        if (~isCZM)
            bPrimeRel = str2num(bTimesPiCZM_s);
            if (isAT1)
                PFtmp.bPrime = bPrimeRel * 3.0 / 8.0;
            else
                PFtmp.bPrime = bPrimeRel * 27.0 / 256.0;
            end
        end
        vals_default = {isHyper, isApproximate, model_s, la_s, l_cD2c_s, df_s, CZM_normalization4AT, bTimesPiCZM_s, CZM_modelName, tauModel};
        [PFtmp, vals_out{i}] = PFtmp.Initialize_Stage1(vals_default);
        PFs{i} = PFtmp;
    end
else
    vec_valsIn = cell(0);
    vec_valsOut = cell(0);
    names = cell(0);
    EPH = '';
    if (isHyper == 1)
        EPH = 'H';
    elseif (isHyper == 0)
        EPH = 'P';
    end

    name_out = [EPH, '_', model_s, '_lcRat_', l_cD2c_s, '_df_', df_s];

    % loading rate
    if ((option >= 30) && (option <= 32))
        if (isHyper == -1)
            isHyper = 1;
        end
        if (isHyper == 1)
            header = '$$ \mathrm{log}_{10}({\bar{\dot{\epsilon}}}_{{\circ}H})  $$';
        else
            header = '$$ \mathrm{log}_{10}({\bar{\dot{\epsilon}}}_{{\circ}P})  $$';
        end
        % start of adding approximate ones
        vec_la = [];
        vec_lapprox = [];
        vec_laNum = [];

        addElliptic = 0;
        if (option > 30)
            drawPlots = 1;
            la_st = -4;
            la_step = 0.25; % use 0.5 for faster generation of data
            la_en = 6; % was 4
            vec_la = la_st:la_step:la_en;
            if (option == 31)
                isApproximate = 0;
            elseif (option == 32)
                isApproximate = 1;
                la_st = 0.5;
                vec_la = la_st:la_step:la_en;
            end
            vals_default = {isHyper, isApproximate, model_s, la_s, l_cD2c_s, df_s, CZM_normalization4AT, bTimesPiCZM_s, CZM_modelName, tauModel};
        else
            isApproximate = 2;
            addElliptic = 1;
            vec_la = [-2, -1, 0, 0.5, 1, 1.5, 2, 2.5];
            la_st_approx = 0.5;
            ai = 0;
            for i = 1:length(vec_la)
                la = vec_la(i);
                if (la >= la_st_approx)
                    ai = ai + 1;
                    vec_lapprox(ai) = la;
                    vec_laNum(ai) = i;
                end
            end
        end
        
        addAPproximate = (isApproximate == 2);
        if (addAPproximate)
            isApproximate = 0;
            vals_default = {isHyper, isApproximate, model_s, la_s, l_cD2c_s, df_s, CZM_normalization4AT, bTimesPiCZM_s, CZM_modelName, tauModel};
        else
            la_st_approx = 10;
        end
        cntr = 0;
        if (CZM_normalization4AT == 0)
            option = option + 5;
        end

        % adding based ones (generally exact ones)
        for i = 1:length(vec_la)
            la = vec_la(i);
            la_sn = num2str(la);
            valsIn = vals_default;
            valsIn{4} = la_sn;
            
            cntr = cntr + 1;
            PFtmp = PF;
            [PFtmp, vals_out{cntr}] = PFtmp.Initialize_Stage1(valsIn);
            names{cntr} = la_sn;
            PFs{cntr} = PFtmp;
            lnums{cntr} = cntr + 1;
            lstys{cntr} = '-';
        end

        for i = 1:length(vec_lapprox)
            la = vec_lapprox(i);
            la_sn = num2str(la);
            valsIn = vals_default;
            valsIn{4} = la_sn;
            valsIn{2} = 1;
            lnum = vec_laNum(i) + 1;
            
            cntr = cntr + 1;
            PFtmp = PF;
            [PFtmp, vals_out{cntr}] = PFtmp.Initialize_Stage1(valsIn);
            names{cntr} = ['none', la_sn, '-', HR_name];
            PFs{cntr} = PFtmp;
            lnums{cntr} = lnum;
            lstys{cntr} = '--';
        end

        % adding elliptic
        if (addElliptic)
            cntr = cntr + 1;
            valsIn = vals_default;
            valsIn{1} = -1;
            PFtmp = PF;
            [PFtmp, vals_out{cntr}] = PFtmp.Initialize_Stage1(valsIn);
            names{cntr} = 'EPF';
            PFs{cntr} = PFtmp;
            lnums{cntr} = 1;
            lstys{cntr} = '-';
        end
        if ((option == 30) || (option == 35))
            lnum_extras{1} = 1;
            lstys_extras{1} = '--';
            names_extras{1} = nameHP_HR{isHyper + 1};
        end
    % damping factor
    elseif ((option >= 40) && (option <= 42))
        header = '$$ f_d $$';
        optionBK = option;
        if (isHyper == -1)
            isHyper = 1;
        end
        % start of adding approximate ones
        if (strcmp(la_s, 'none') == 1)
            la_s = '1';
        end
        la = num2str(la_s);
        if (length(la) == 0)
            return;
        end
        if (la <= -100)
            return;
        end
        ldfApprox = [];
        dfApprox = [];
        vec_dfNum = [];

        addElliptic = 0;
        addParabolic = 0;
        

        if (optionBK > 40)
            drawPlots = 0;
            ldf = -3:0.25:2;
%            ldf = -3:1:2;
%            drawPlots = 1;
            df = power(10, ldf);
            if (optionBK == 41)
                isApproximate = 0;
            elseif (optionBK == 42)
                isApproximate = 1;
                if ((la < 0.5) || (isHyper == 1))
                    return;
                end
            end
            vals_default = {isHyper, isApproximate, model_s, la_s, l_cD2c_s, df_s, CZM_normalization4AT, bTimesPiCZM_s, CZM_modelName, tauModel};
        else
            addApprox = (isApproximate ~= 0);
            if (la <= 0)
                addApprox = 0;
            end
            if (isHyper == 1)
                addElliptic = 1;
                addParabolic = 1;
            elseif (isHyper == 0)
                addElliptic = 1;
            end
%            ldf = -2.0:.5:1.5;
            ldf = [-2, -1, 0, 1, 2]; 
            df = power(10, ldf);
            
            if (addApprox)
                if (isHyper == 0)
                    ldfApprox = ldf;
                    dfApprox = df;
                    vec_dfNum = 1:length(df);
                else
                    sz_df = length(df);
                    ldfApprox = ldf(sz_df);
                    dfApprox = df(sz_df);
                    vec_dfNum = -1;
                end
            end
        end
        if (isApproximate == 2)
            isApproximate = 0;
            vals_default = {isHyper, isApproximate, model_s, la_s, l_cD2c_s, df_s, CZM_normalization4AT, bTimesPiCZM_s, CZM_modelName, tauModel};
        end
        if (CZM_normalization4AT == 0)
            option = option + 5;
        end
        cntr = 0;

        % adding based ones (generally exact ones)
        sz = length(df);
        sz_approx = length(dfApprox);

        lnumBase = [];
        for i = 1:sz
            dfv = df(i);
            df_s = num2str(dfv);
            ldfv = ldf(i);
            valsIn = vals_default;
            valsIn{6} = df_s;
            cntr = cntr + 1;
            PFtmp = PF;
            [PFtmp, vals_out{cntr}] = PFtmp.Initialize_Stage1(valsIn);
            names{cntr} = df_s;
            PFs{cntr} = PFtmp;
            lnums{cntr} = cntr + 1;
            lnumBase(i) = lnums{cntr};
            lstys{cntr} = '-';
        end

        for i = 1:sz_approx
            dfv = dfApprox(i);
            df_s = num2str(dfv);
            ldfv = ldfApprox(i);

            valsIn = vals_default;
            valsIn{6} = df_s;
            valsIn{2} = 1;
            lnum = vec_dfNum(i) + 1;
            nm = ['none', df_s, '-', HR_name];
            if (lnum <= 0)
                lnum = 1;
                nm = nameHP_HR{2};
            end
            cntr = cntr + 1;
            PFtmp = PF;
            [PFtmp, vals_out{cntr}] = PFtmp.Initialize_Stage1(valsIn);
            names{cntr} = nm;
            PFs{cntr} = PFtmp;
            lnums{cntr} = lnum;
            lstys{cntr} = '--';
        end

        % adding elliptic
        if (addElliptic)
            cntr = cntr + 1;
            valsIn = vals_default;
            valsIn{1} = -1;
            PFtmp = PF;
            [PFtmp, vals_out{cntr}] = PFtmp.Initialize_Stage1(valsIn);
            names{cntr} = 'EPF';
            PFs{cntr} = PFtmp;
            lnums{cntr} = 1;
            lstys{cntr} = '-';
        end

        % adding parabolic
        if (addParabolic)
            for i = 1:sz
                dfv = df(i);
                df_s = num2str(dfv);
                ldfv = ldf(i);
                valsIn = vals_default;
                valsIn{1} = 0;
                valsIn{6} = df_s;
                cntr = cntr + 1;
                PFtmp = PF;
                [PFtmp, vals_out{cntr}] = PFtmp.Initialize_Stage1(valsIn);
                names{cntr} = ['nonePPF-',df_s];
                PFs{cntr} = PFtmp;
                lnums{cntr} = lnumBase(i);
                lstys{cntr} = ':';
            end
        end
        if (optionBK == 40)
            if (isHyper == 0) % PPF
                lnum_extras{1} = 101;
                lstys_extras{1} = '--';
                names_extras{1} = nameHP_HR{isHyper + 1};
            else
                lnum_extras{1} = 101;
                lstys_extras{1} = ':';
                names_extras{1} = nameHP{1};
            end
        end
        
    % wave speed ratio
    elseif ((option >= 50) && (option <= 52))
        header = '$$ k  $$';
        optionBK = option;
        if (isHyper == -1)
            isHyper = 1;
        end
        % start of adding approximate ones
        if (strcmp(la_s, 'none') == 1)
            la_s = '1';
        end
        la = num2str(la_s);
        if (length(la) == 0)
            return;
        end
        if (la <= -100)
            return;
        end
        lcRatApprox = [];
        cRatApprox = [];
        vec_cRatNum = [];

        addElliptic = 0;
        addParabolic = 0;

        if (optionBK > 50)
            drawPlots = 0;
            lcRat = -3:0.25:2;
            cRat = power(10, lcRat);
            if (optionBK == 51)
                isApproximate = 0;
            elseif (optionBK == 52)
                isApproximate = 1;
                if (la < 0.5)
                    return;
                end
            end
            vals_default = {isHyper, isApproximate, model_s, la_s, l_cD2c_s, df_s, CZM_normalization4AT, bTimesPiCZM_s, CZM_modelName, tauModel};
        else
            addApprox = (isApproximate ~= 0);
            if (la <= 0)
                addApprox = 0;
            end
            if (additionalLines4_50 == 0)
                addApprox = 0;
            end
            if (isHyper == 1)
                addElliptic = 1;
                addParabolic = 1;
                if (additionalLines4_50 == 0)
                    addParabolic = 0;
                end
            elseif (isHyper == 0)
                addElliptic = 1;
            end
%            ldf = -2.0:.5:1.5;
            lcRat = [-2, -1, 0, 1, 2]; 
            if (isHyper == 1)
                lcRat = [-1, 0, 1, 2]; 
            end
            cRat = power(10, lcRat);
            
            if (addApprox)
                if (1) %(isHyper == 0)
                    lcRatApprox = lcRat;
                    cRatApprox = cRat;
                    vec_dfNum = 1:length(cRat);
                else
                    sz_cRat = length(cRat);
                    ldfApprox = ldf(sz_cRat);
                    dfApprox = df(sz_cRat);
                    vec_cRatNum = -1;
                end
            end
            if ((optionBK == 50) && (addParabolic == 1))
                if (isHyper == 0) % PPF
                    lnum_extras{1} = 101;
                    lstys_extras{1} = '--';
                    names_extras{1} = nameHP_HR{isHyper + 1};
                else
                    lnum_extras{1} = 101;
                    lstys_extras{1} = ':';
                    names_extras{1} = nameHP{1};
                end
            end
        end
        if (isApproximate == 2)
            isApproximate = 0;
            vals_default = {isHyper, isApproximate, model_s, la_s, l_cD2c_s, df_s, CZM_normalization4AT, bTimesPiCZM_s, CZM_modelName, tauModel};
        end
        if (CZM_normalization4AT == 0)
            option = option + 5;
        end
        cntr = 0;

        % adding based ones (generally exact ones)
        sz = length(cRat);
        sz_approx = length(cRatApprox);

        lnumBase = [];
        for i = 1:sz
            cRatv = cRat(i);
            cRat_s = num2str(cRatv);
            lcRatv = lcRat(i);
            clRat_s = num2str(lcRatv);
            valsIn = vals_default;
            valsIn{5} = clRat_s;
            cntr = cntr + 1;
            PFtmp = PF;
            [PFtmp, vals_out{cntr}] = PFtmp.Initialize_Stage1(valsIn);
            names{cntr} = cRat_s;
            PFs{cntr} = PFtmp;
            lnums{cntr} = cntr + 1;
            lnumBase(i) = lnums{cntr};
            lstys{cntr} = '-';
        end

        for i = 1:sz_approx
            cRatv = cRatApprox(i);
            cRat_s = num2str(cRatv);
            lcRatv = lcRatApprox(i);
            clRat_s = num2str(lcRatv);
            valsIn = vals_default;
            valsIn{5} = clRat_s;
            valsIn{2} = 1;
            lnum = vec_dfNum(i) + 1;
            nm = ['none', df_s, '-', HR_name];
%            if (lnum <= 0)
%                lnum = 1;
%                nm = ['HPF-high-rate'];
%            end
            cntr = cntr + 1;
            PFtmp = PF;
            [PFtmp, vals_out{cntr}] = PFtmp.Initialize_Stage1(valsIn);
            names{cntr} = nm;
            PFs{cntr} = PFtmp;
            lnums{cntr} = lnum;
            lstys{cntr} = '--';
        end

        % adding elliptic
        if (addElliptic)
            cntr = cntr + 1;
            valsIn = vals_default;
            valsIn{1} = -1;
            PFtmp = PF;
            [PFtmp, vals_out{cntr}] = PFtmp.Initialize_Stage1(valsIn);
            names{cntr} = 'EPF';
            PFs{cntr} = PFtmp;
            lnums{cntr} = 1;
            lstys{cntr} = '-';
        end

        % adding parabolic
        if (addParabolic)
            for i = 1:sz
                cRatv = cRat(i);
                cRat_s = num2str(cRatv);
                lcRatv = lcRat(i);
                clRat_s = num2str(lcRatv);
                valsIn = vals_default;
                valsIn{5} = clRat_s;
                valsIn{1} = 0;
                cntr = cntr + 1;
                PFtmp = PF;
                [PFtmp, vals_out{cntr}] = PFtmp.Initialize_Stage1(valsIn);
                names{cntr} = ['nonePPF-',df_s];
                PFs{cntr} = PFtmp;
                lnums{cntr} = lnumBase(i);
                lstys{cntr} = ':';
            end
        end
        if ((optionBK == 50) && (addParabolic == 1))
            if (isHyper == 0) % PPF
                lnum_extras{1} = 101;
                lstys_extras{1} = '--';
                names_extras{1} = nameHP_HR{isHyper + 1};
            else
                lnum_extras{1} = 101;
                lstys_extras{1} = '--';
                names_extras{1} = nameHP_HR{2}; %nameHP{2};
                lnum_extras{2} = 101;
                lstys_extras{2} = ':';
                names_extras{2} = nameHP{1};
            end
        end
    end
end

num = length(PFs);
for i = 1:num
    PFs{i} = PFs{i}.Compute();
end
vecs_names = PFs{1}.vecs_names;
scalar_names = PFs{1}.scalar_names;
num_vec = length(vecs_names);
num_scalar = length(scalar_names);

folderName = ['option_', num2str(option)];
if (length(name_out) > 0)
    folderName = [folderName, '_', name_out];
else
    folderName = [folderName, '_ap_', num2str(ap), '_isHyper_', num2str(isHyper)];
end
[status,message,messageid] = mkdir(folderName);
fnws = folderName;
folderName = [folderName, '/'];

if ((useColor4Lines == 0) && (num > 4))
    useColor4Lines = 1;
end
linePropSpecified = (length(lnums) > 0);

colorNum = 1;
isModel = 0;
includeBlack = 0;
if (option >= 30)
    includeBlack = 1;
end
lc = getColors(colorNum, isModel, includeBlack);



lthCol = 2.5;
ls = {'-', '--','-.', ':'};
lthk = [2, 2, 2, lthCol];


has_col_out = (~isempty(vals_out));
num_col = length(vals_default);

fn = [folderName, 'scalars.csv'];
fid = fopen(fn, 'w');
fprintf(fid, 'runName');
if (has_col_out)
    for i = 1:num_col
        fprintf(fid, ',%s', name_vals{i});
    end
end

for i = 1:num_scalar
    fprintf(fid, ',%s', scalar_names{i});
end
fprintf(fid, '\n');
for ri = 1:num
    nm = names{ri};
    if contains(nm, 'none')
        nm = strrep(nm, 'none', '');    
    end
    fprintf(fid, '%s', nm);
    if (has_col_out)
        for i = 1:num_col
            fprintf(fid, ',%s', vals_out{ri}{i});
        end
    end
    for i = 1:num_scalar
        fprintf(fid, ',%g', PFs{ri}.scalars(i));
    end
    fprintf(fid, '\n');
end
fclose(fid);

if (drawPlots)
    printfig = ~(((option == 31) || (option == 32) || (option == 36) || (option == 37)) && (isHyper == 1));
    % forming names and clrs for the figure
    % actual goes to the figure
    actual_clrs = cell(0);
    actual_th = cell(0);
    actual_sty = cell(0);
    % leg are those that appear in the legend
    leg_names{1} = header;
    leg_clrs{1} = 'none';
    leg_th{1} = 1;
    leg_sty{1} = 'none';
    leg_cntr = 1;
    
    for ri = 1:num
        legname = names{ri};
        add2Leg = (~contains(legname, 'none'));

        lcv = 'k';
        lsv = '-';
        lthv = 2;

        if (linePropSpecified)
            lcv = lc{lnums{ri}};
            lsv = lstys{ri};
        else
            if (useColor4Lines)
                lcv = lc{ri};
            else
                lsv = ls{ri};
%                lthv = lth{ri};
            end
        end
        if (lsv == ':')
            lthv = lthCol;
        end

        actual_clrs{ri} = lcv;
        actual_th{ri} = lthv;
        actual_sty{ri} = lsv;
        
        if (add2Leg)
            leg_cntr = leg_cntr + 1;
            leg_names{leg_cntr} = legname;
            leg_clrs{leg_cntr} = lcv;
            leg_th{leg_cntr} = lthv;
            leg_sty{leg_cntr} = lsv;
        end
    end

    numExtras = length(lnum_extras);
    for ri = 1:numExtras
        legname = names_extras{ri};
%        add2Leg = (~contains(legname, 'none'));

        lcv = lc{lnum_extras{ri}};
        lsv = lstys_extras{ri};
        lthv = 2;
        if (lsv == ':')
            lthv = lthCol;
        end
        leg_cntr = leg_cntr + 1;
        leg_names{leg_cntr} = legname;
        leg_clrs{leg_cntr} = lcv;
        leg_th{leg_cntr} = lthv;
        leg_sty{leg_cntr} = lsv;
    end
    
    lfs = 13;
    if (leg_cntr > 30)
        lfs = 6;
    elseif (leg_cntr > 20)
        lfs = 7;
    elseif (leg_cntr > 10)
        lfs = 10;
    end        

    PFtm = PF;
    vecs_names_latex = PFtm.vecs_names_latex;

    for vi = 1:num_vec - 1
        vec_name = vecs_names{vi};
        fnbase = [folderName, fnws, '_eps_vs_', vec_name];
        fg = figure(vi);
        set(fg,'defaultLegendAutoUpdate', 'off');

        for lri = 1:leg_cntr
            plot([nan], [nan], 'Color', leg_clrs{lri}, 'LineStyle', leg_sty{lri}, 'LineWidth', leg_th{lri});
            hold on;
        end
        lg = legend(leg_names, 'FontSize', lfs, 'Interpreter', 'latex');
        legend('boxoff');


        for ri = 1:num
            lcv = actual_clrs{ri};
            lthv = actual_th{ri};
            lsv = actual_sty{ri};
            
            x = PFs{ri}.epsilon_p;
            y = PFs{ri}.y_vecs{vi};
            plot(x, y, 'Color', lcv, 'LineStyle', lsv, 'LineWidth', lthv);
            hold on;
        end
    
        xh = get(gca, 'XLabel');
        set(xh, 'String', '$$ \bar{\epsilon} $$', 'FontSize', labsz, 'VerticalAlignment','Top', 'Interpreter', 'latex');
        yh = get(gca, 'YLabel');
        set(yh, 'String', vecs_names_latex{vi}, 'FontSize', labsz, 'VerticalAlignment','Bottom', 'Interpreter', 'latex');
    
        print('-dpng', [fnbase, '.png']);
        if (printfig)
            savefig([fnbase, '.fig']);
        end
    end
    close('all');
    if (plot_D_xAxis)
        for vi = 1:num_vec
            if (vi == 2)
                continue;
            end
            vec_name = vecs_names{vi};
            fnbase = [folderName, fnws, '_D_vs_', vec_name];
            fg = figure(vi + 10);
            set(fg,'defaultLegendAutoUpdate', 'off');
    
            for lri = 1:leg_cntr
                plot([nan], [nan], 'Color', leg_clrs{lri}, 'LineStyle', leg_sty{lri}, 'LineWidth', leg_th{lri});
                hold on;
            end
            lg = legend(leg_names, 'FontSize', lfs, 'Interpreter', 'latex');
            legend('boxoff');
    
    
            for ri = 1:num
                lcv = actual_clrs{ri};
                lthv = actual_th{ri};
                lsv = actual_sty{ri};
                
                x = PFs{ri}.D_vec;
                y = PFs{ri}.y_vecs{vi};
                plot(x, y, 'Color', lcv, 'LineStyle', lsv, 'LineWidth', lthv);
                hold on;
            end
        
            xh = get(gca, 'XLabel');
            set(xh, 'String', '$$ D $$', 'FontSize', labsz, 'VerticalAlignment','Top', 'Interpreter', 'latex');
            yh = get(gca, 'YLabel');
            set(yh, 'String', 'ylabel', 'FontSize', labsz, 'VerticalAlignment','Bottom', 'Interpreter', 'latex');
        
            print('-dpng', [fnbase, '.png']);
            if (printfig)
                savefig([fnbase, '.fig']);
            end
        end
    end
    close('all');


end
fclose('all');
