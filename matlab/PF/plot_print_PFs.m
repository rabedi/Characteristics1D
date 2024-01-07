function plot_print_PFs(option, ap, isHyper)
labsz = 25; % x, y label font size
%option 0 to 10 are for comparing various PF formulations
%option 0: AT1
%option 1: AT1, AT2
%option 2: CZM-Wu, CZM-Lorentz
%option 5: AT1, AT2, CZM-Wu, CZM-Lorentz
%option 6: AT1, CZM-Wu, CZM-Lorentz


%option 10: CZM-Wu: for all TSR models (linear, ...)

%option 20: CZM-Wu: linear model, for different b'
%%%%%%%%%%%%%%%%%%%
% ap is the loading rate:
% if > 0, it's real, if < 0, for H and P PF, it computes the asymptotic
% high loading rate solution

% isHyper
%       -1: Elliptic
%        0: Parabolic
%        1: Hyperbolic

if nargin < 1
    option = 5; %5, 2, 10, 20
end

if nargin < 2
    ap = -1; %-1; %100;
end


if nargin < 3
    isHyper = -1;
end

% each curve is distniguished by a color
useColor4Lines = 1;
names = cell(0);
PFs = cell(0);
if (option < 10)
    CZM_normalization4AT1_2 = 1;
    xis = [];
    omegaCZMs = [];
    if (option == 0) % AT1
        xis = [1];
        omegaCZMs = [0];
        names = {'AT1'};
        CZM_normalization4AT1_2 = 0;
    elseif (option == 1) T
        xis = [1, 0];
        omegaCZMs = [0, 0];
        names = {'AT1', 'AT2'};
        CZM_normalization4AT1_2 = 0;
    elseif (option == 2) % CZM-Wu, CZM-Lorentz
        xis = [2, 1];
        omegaCZMs = [1, 1];
        names = {'CZM-Wu', 'CZM-Lorentz'};
    elseif (option == 5) % AT1, AT2, CZM-Wu, CZM-Lorentz
        xis = [1, 0, 2, 1];
        omegaCZMs = [0, 0, 1, 1];
        names = {'AT1', 'AT2', 'CZM-Wu', 'CZM-Lorentz'};
    elseif (option == 6) % AT1, CZM-Wu, CZM-Lorentz
        xis = [1, 2, 1];
        omegaCZMs = [0, 1, 1];
        names = {'AT1', 'CZM-Wu', 'CZM-Lorentz'};
    end
    num = length(xis);
    PFs = cell(num, 1);
    for i = 1:num
        PFtmp = PF;
        PFtmp.ap = ap;
        PFtmp.isHyper = isHyper;
        PFtmp.CZM_normalization4AT1_2 = CZM_normalization4AT1_2;
        PFtmp.xi = xis(i);
        PFtmp.omegaCZM = omegaCZMs(i);
        PFtmp.CZM_model_name = 'Linear';
        PFtmp.bPrime = -1;
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
elseif (option == 20)
    names = {'2', '1', '0.5'}; % bPrime * pi
    num = length(names);
    for i = 1:num
        PFtmp = PF;
        PFtmp.ap = ap;
        PFtmp.isHyper = isHyper;
        PFtmp.CZM_normalization4AT1_2 = 1;
        PFtmp.xi = 2;
        PFtmp.omegaCZM = 1;
        PFtmp.CZM_model_name = 'Linear';
        PFtmp.bPrime = str2num(names{i}) / pi;
        PFs{i} = PFtmp;
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

forlderName = ['option_', num2str(option), '_ap_', num2str(ap), '_isHyper_', num2str(isHyper)];
[status,message,messageid] = mkdir(forlderName);
forlderName = [forlderName, '/'];

if ((useColor4Lines == 0) && (num > 4))
    useColor4Lines = 1;
end

lc = cell(0);
cntr = 1;
lc{cntr} = [0 0 0]; % black
cntr = cntr + 1;
lc{cntr} = [1	0	0]; %red'
cntr = cntr + 1;
lc{cntr} = [1	102/255	0]; %orange'
cntr = cntr + 1;
%lc{cntr} = [0.5	0	1]; % purple 
%cntr = cntr + 1;
lc{cntr} = [0	0	1]; % blue 'b';
cntr = cntr + 1;
lc{cntr} = [0 135/255 0]; %green
cntr = cntr + 1;
lc{cntr} = [0	1	1]; % teal [0	0	0.5]; % dark blue
cntr = cntr + 1;
lc{cntr} = [0.5	0.25	0]; % brown
cntr = cntr + 1;
lc{cntr} = [0.5	0.5	0.5]; %'dark_gray2'
cntr = cntr + 1;
lc{cntr} = [0.75	0.75	0.75]; %'gray2'
cntr = cntr + 1;
lc{cntr} = [1	1	0]; %yellow [0.5	0.5	0]; % olive , [0	0.5	0.25]; %green blue
cntr = cntr + 1;

ls = {'-', '--','-.', ':'};
lthk = [2, 2, 2, 3];


fn = [forlderName, 'scalars.csv'];
fid = fopen(fn, 'w');
fprintf(fid, 'runName');
for i = 1:num_scalar
    fprintf(fid, ',%s', scalar_names{i});
end
fprintf(fid, '\n');
for ri = 1:num
    fprintf(fid, '%s', names{ri});
    for i = 1:num_scalar
        fprintf(fid, ',%g', PFs{ri}.scalars(i));
    end
    fprintf(fid, '\n');
end
fclose(fid);

for vi = 1:num_vec
    vec_name = vecs_names{vi};
    fnbase = [forlderName, vec_name];
    figure(vi);
    for ri = 1:num
        lcv = 'k';
        lsv = '-';
        lthv = 2;
        if (useColor4Lines)
            lcv = lc{ri};
        else
            lsv = ls{ri};
            lthv = lth{ri};
        end
        x = PFs{ri}.epsilon_p;
        y = PFs{ri}.y_vecs{vi};
        plot(x, y, 'Color', lcv, 'LineStyle', lsv, 'LineWidth', lthv);
        hold on;
    end
    lg = legend(names, 'FontSize', 12);
    legend('boxoff');

    xh = get(gca, 'XLabel');
    set(xh, 'String', 'xlabel', 'FontSize', labsz, 'VerticalAlignment','Top', 'Interpreter', 'latex');
    yh = get(gca, 'YLabel');
    set(yh, 'String', 'ylabel', 'FontSize', labsz, 'VerticalAlignment','Bottom', 'Interpreter', 'latex');

    print('-dpng', [fnbase, '.png']);
    savefig([fnbase, '.fig']);
end
fclose('all');