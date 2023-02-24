function plot_TimeSequences_subdomains(dim, ringProb)
    root = '../run_0';
if nargin < 1
    dim = 1;
end
if nargin < 2
    ringProb = 0;
end

if (strcmp(root, '') == 0)
    rt = [root, '/'];
else
    rt = '';
end


labels = {'time',	'sigmaL', 'sigmaR', 'vL',	'vR',  'uL', 'uR', 'delu', ...
    'pL', 'pR', 'eL', 'del_e', 'diss', ...,
    'damage', 'damageSource', 'maxEffdelu', 'v_r', 'v_r_source', 'sigma_theta_source'};
ln = length(labels);
if (ringProb == 0)
    ln = ln - 3;
end

pos = 0;
fileName{1} = [rt, 'xPos_dl_0_dR_0_InterfaceRawFinalSln.txt'];
lg{1} = '0-0';
fileName{2} = [rt, 'xPos_dl_0_dR_1_InterfaceRawFinalSln.txt'];
lg{2} = '0-1';
fileName{3} = [rt, 'xPos_dl_1_dR_2_InterfaceRawFinalSln.txt'];
lg{3} = '1-2';
fileName{4} = [rt, 'xPos_dl_2_dR_2_InterfaceRawFinalSln.txt'];
lg{4} = '2-2';

sz = 0;
for ii = 1:4
    fn = fileName{ii};
    fid = fopen(fn, 'r');
    if (fid < 0)
        continue;
    end
    sz = sz + 1;
    fclose(fid);
    [data{sz}, delu{sz}, delv{sz}, sigmaL{sz}, sigmaR{sz}, diss{sz}, del_e{sz}, exists] = loadTimeSequence(fn, dim, ringProb);
    lgs{sz} = lg{ii};
end

figure(1);
clf
for i = 1:sz
   du = delu{i};
   sg = sigmaR{i};
   plot(du, sg);
   leg{i} = num2str(i - 1);
   hold on;
end
lgh = legend(lgs, 'FontSize', 9);
legend('boxoff');
plotName = [rt, 'plot_delu_sigma'];
print('-dpng', [plotName, '.png']);
savefig([plotName, '.fig']);
clf

[a, ln] = size(data{1});
for j = 2:ln
    figure(j);
    for k = 1:sz
        dataOut = data{k};
        times = dataOut{1};
        dat = dataOut{j};
        y = dat(:, 1);
        plot(times, y);
        hold on;
    end
    lgh = legend(lgs, 'FontSize', 9);
    legend('boxoff');
    plotName = [rt, 'plot', labels{j}];
    print('-dpng', [plotName, '.png']);
    savefig([plotName, '.fig']);
end

fclose('all');
close('all');