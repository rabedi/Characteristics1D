function compare_space_spacetime_results(isSummary, serialNums, versionNums, timeSerialNumbers, timeStamp, subdomainNo)
if nargin < 1
    isSummary = 1;
end
if nargin < 2
    serialNums = [0, 1];
end
if nargin < 3
    versionNums = [11, 12, 13];
end
if nargin < 4
    timeSerialNumbers = [];
end
if nargin < 5
    timeStamp = '';
end
if nargin < 6
    subdomainNo = 0;
end
% if < 0, does not print space files
if nargin < 4
    maxTimeStep = 100000000;
end
if nargin < 5
    step = -1;
end
sd = ['sd_', num2str(subdomainNo), '_'];

sz_serialNums = length(serialNums);
sz_versionNums = length(versionNums);


cntr = 1;
for sii = 1:sz_serialNums
    si = serialNums(sii);
    sis = num2str(si);
    for vii = 1:sz_versionNums
        vi = versionNums(vii);
        leg{cntr} = '';
        if (vi >= 0)
            vis = num2str(vi);
            leg{cntr} = ['i=', sis, ',v=',vis];
        end
        cntr = cntr + 1;
    end
end

dir4Plots = '../_plots';
[status,message,messageid] = mkdir(dir4Plots);
dir4Plots = [dir4Plots, '/'];

lfs = 18;
cntr = 0;
last = sz_serialNums * sz_versionNums;
for sii = 1:sz_serialNums
    si = serialNums(sii);
    sis = num2str(si);
    for vii = 1:sz_versionNums
        vi = versionNums(vii);
        if (vi >= 0)
            vis = num2str(vi);
            vis = ['_V_', vis];
        else
            vis = '';
        end
        cntr = cntr + 1;
        rt = ['../run', vis, '_', sis, '/'];
        if (isSummary)
            prename = [sd, '_Summary'];
            startCol = 2;
        else
            prename = [sd, 'BulkInterface_tPos_', num2str(timeSerialNumbers(sii, vii))];
            startCol = 3;
        end
        fn = [rt , prename, '.txt'];
        fid = fopen(fn, 'r');
        if (fid > 0)
            dwh = datawheader;
            dwh = dwh.read(fid);
            fclose(fid);

            x = dwh.dataMat(:,startCol);
            for j = startCol + 1:dwh.nFieldCols
                figure(j);
                y = dwh.dataMat(:, j);
                plot(x, y);
                hold on;
                if (cntr == last)
                    xh = get(gca, 'XLabel');
                    set(xh, 'String', ['$$ ', dwh.fieldLatexNames{startCol}, ' $$'], 'FontSize', lfs, 'VerticalAlignment','Top', 'Interpreter', 'latex');
    
                    yh = get(gca, 'YLabel');
                    set(yh, 'String', ['$$ ', dwh.fieldLatexNames{j}, ' $$'], 'FontSize', lfs, 'VerticalAlignment','Bottom', 'Interpreter', 'latex');

                    if (isSummary)
                        name = [dir4Plots, prename, 'plot_', num2str(startCol, '%02d'), '_', num2str(j, '%02d'), '_', dwh.fieldNames{j}];
                    else
                        name = [dir4Plots, prename, 'plot_', '_ti_', num2str(timeSerialNumbers(1, 1)), '_', num2str(startCol, '%02d'), '_', num2str(j, '%02d'), '_', dwh.fieldNames{j}];
                        title(['time stamp = ', num2str(timeStamp)]);
                    end
                    legend(leg, 'FontSize', 10);
                    legend('boxoff');
                    savefig([name, '.fig']);
                    print('-dpng', [name, '.png']);
                end
            end
        end
    end
end
fclose('all');
close('all');