function Compare_OneFld_Stats(fld, ddStr, doWeibull, doGauss)
lfs = 14;
labsz = 25;
if (nargin < 1)
    fld = 1;% 2
end
if (nargin < 2)
    ddStr = '0.5';
end
if (nargin < 3)
    doWeibull = -1;
end
if (nargin < 4)
    doGauss = -1;
end
if (doGauss == -1)
    doGauss = (fld == 1);
end
if (doWeibull == -1)
    doWeibull = (fld == 2);
end

fldcstr = ['fld', num2str(fld)];

shapeStrs = {'1', '1.5', '2', '3', '4'};
n_shape = length(shapeStrs);
llcStrs = {'-4.5', '-4', '-3.5', '-3', '-2.5', '-2', '-1.5', '-1', '-0.5'};
%llcStrs = {'-4', '-3.5', '-3', '-2.5', '-2', '-1.5', '-1', '-0.5'};

n_llcStr = length(llcStrs);
for llci = 1:n_llcStr
    llcs(llci) = str2num(llcStrs{llci});
end

for si = 1:n_shape
    shapeStr = shapeStrs{si};
    for llci = 1:n_llcStr
        llcStr = llcStrs{llci};
        fileName = ['stat/', fldcstr, '/shape', shapeStr, '_dd', ddStr, '_llc', llcStr, '_', fldcstr, '_stat.txt'];
        [vals, names, n] = ReadStatFile(fileName);
        for i = 1:n
            data{i}{si}(llci) = vals(i);
        end
    end
end

lc{1} = 'k';
lc{2} = 'r';
lc{3} = 'b';
lc{4} = 'c';
lc{5} = 'g';


[status,msg,msgID] = mkdir('plots');
base0 = ['plots/', fldcstr];
[status,msg,msgID] = mkdir(base0);
for i = 1:n
    figure(1);
    clf
    name = ['t', num2str(i), names{i}];
    base1 = [base0, '/', name];
    [status,msg,msgID] = mkdir(base1);
    fnbase = [base1, '/plot_', fldcstr, '_', name];
    for si = 1:n_shape
        y = data{i}{si};
        plot(llcs, y, 'Color', lc{si}, 'LineWidth', 2);
        hold on;
    end
    lg = legend(shapeStrs, 'FontSize', lfs, 'Interpreter', 'latex');
    legend('boxoff');
    xh = get(gca, 'XLabel');
    set(xh, 'String', '$$ \mathrm{log}_{10}(l_{\mathrm{cor}}) $$', 'FontSize', labsz, 'VerticalAlignment','Top', 'Interpreter', 'latex');
    yh = get(gca, 'YLabel');
    set(yh, 'String', 'ylabel', 'FontSize', labsz, 'VerticalAlignment','Bottom', 'Interpreter', 'latex');

    print('-dpng', [fnbase, '.png']);
    savefig([fnbase, '.fig']);
end