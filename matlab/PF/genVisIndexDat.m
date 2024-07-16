classdef genVisIndexDat
    properties
        % plots columns versus one x axis
        colNamesLatex = cell(0);
        xLabelLatex = '$$ \mathrm{xlabel} $$';

        
        % 5 groups are supported
        % indexed as i0, i1, i2, i3, i4, i5
        % data{i0}{i1}{i2}{i3}{i4}{i5} -> [nRow x ncol] data
        % 
        data;
        % again indexed as {i0} .... stored x values, size is nRow
        xVals;
        
        % 
        % size 5 cell, corresponding to groups, for each group, the
        % available names are given in latex
        groupNames = cell(0);
       
        lfs = 12;
        labsz = 25; % x, y label font size
    end
    methods
        % from 5 groups:
        % clrPos -> which group goes to color

        % lineStylePos -> which group goes to line style

        % which group goes to markers

        function Visualize(obj, optionNo, addedFolderFileName, starts, ends, clrPos, lineStylePos, markerPos, added_indices)
            if (nargin < 9)
                added_indices = cell(0);
            end
            option121_hack = ((optionNo == 121) || (optionNo == 621));
            option200_df_hack = ((optionNo == 200) || (optionNo == 700));
            optioNo_str = num2str(optionNo);
            xLabel = obj.xLabelLatex;

            markers = {'+', 'o', 's', 'p', '*'};
    
            inds = cell(0);
            i = 1;
            for i1 = starts(1):ends(1)
                for i2 = starts(2):ends(2)
                    for i3 = starts(3):ends(3)
                        for i4 = starts(4):ends(4)
                            for i5 = starts(5):ends(5)
                                inds{i} = [i1, i2, i3, i4, i5];
                                i = i + 1;
                            end
                        end
                    end
                end
            end
            numPtsBeforeAdded = length(inds);
            sz2 = length(added_indices);
            for j = 1:sz2
                indv = added_indices{j};
                i1 = indv(1);
                i2 = indv(2);
                i3 = indv(3);
                i4 = indv(4);
                i5 = indv(5);
                inds{i} = [i1, i2, i3, i4, i5];
                i = i + 1;
            end
            numPts = length(inds);
    
            lthCol = 2.5;
            ls = {'-', '--','-.', ':'};
            lthg = 2;
            lthk = [lthg, lthg, lthg, lthCol];
    
            nolc = [0.5	0.5	0.5]; %'dark_gray2'
            clrwh = [1, 1, 1];
            
            colorNum = 0;
            isModel = 0;
            includeBlack = ~option200_df_hack;
            lc = getColors(colorNum, isModel, includeBlack);
            if 0
                lc = cell(0);
                cntr = 0;
                if (~option200_df_hack)
                    cntr = cntr + 1;
                    lc{cntr} = [0 0 0]; % black
                end
                cntr = cntr + 1;
                lc{cntr} = [1	0	0]; %red'
                cntr = cntr + 1;
                lc{cntr} = [0	0	1]; % blue 'b';
                cntr = cntr + 1;
                lc{cntr} = [0 135/255 0]; %green
                cntr = cntr + 1;
                lc{cntr} = [1	102/255	0]; %orange'
                cntr = cntr + 1;
                lc{cntr} = [0	1	1]; % teal [0	0	0.5]; % dark blue
                cntr = cntr + 1;
                lc{cntr} = [0.5	0.25	0]; % brown
                cntr = cntr + 1;
                lc{cntr} = 1/255*[255	0	255]; % magenta
                cntr = cntr + 1;
                lc{cntr} = 1/255*[255	128	192]; % 'rosy_pink'
                cntr = cntr + 1;
                lc{cntr} = [0.5	0.5	0.5]; %'dark_gray2'
                cntr = cntr + 1;
                lc{cntr} = [0	0.5	0.25]; %green blue
                cntr = cntr + 1;
                lc{cntr} = [0.75	0.75	0.75]; %'gray2'
                cntr = cntr + 1;
                lc{cntr} = [1	1	0]; %yellow , [0	0.5	0.25]; %green blue
                cntr = cntr + 1;
                lc{cntr} = [0.5	0	1]; % purple 
                cntr = cntr + 1;
                lc{cntr} = 1/255 * [203 0   51]; % red2
                cntr = cntr + 1;
                lc{cntr} = [0.5	0.5	0]; % olive 
                cntr = cntr + 1;
                lc{cntr} = 1/255 * [64	128	128]; % blue2
                cntr = cntr + 1;
                lc{cntr} = 1/255 * [255	128	192]; % rosy_pink
                cntr = cntr + 1;
                lc{cntr} = 1/255 * [255	128	128]; % peach
                cntr = cntr + 1;
                lc{cntr} = 1/255 * [128	0	64]; % arghavani
                cntr = cntr + 1;
                lc{cntr} = 1/255 * [128	0	128]; % purple2
            end    
    
            optionNo_str = num2str(optionNo);
    
            folderName = ['plot_convRate_opt_', optioNo_str];
            [status,message,messageid] = mkdir(folderName);
            hasAddedName = ((~isempty(addedFolderFileName)) && (strcmp(addedFolderFileName, 'none') == 0));
            fnBase0 = folderName;
            if (hasAddedName)
                fnBase0 = [fnBase0, '_', addedFolderFileName];
                folderName = [folderName, '/', addedFolderFileName];
                [status,message,messageid] = mkdir(folderName);
            end
    
             ncol = length(obj.colNamesLatex);

            % plotting
            for col = 1:ncol
                icol = num2str(col);
                yLabel = obj.colNamesLatex{col};
                fg = figure(col);
                set(fg,'defaultLegendAutoUpdate', 'off');
            
                % creating legend
                x = [nan];
                y = [nan];
                cntr = 1;
                plot(x, y, 'Color', 'none');
                hold on;
                leg{cntr} = 'header';
            
                if (clrPos > 0)
                    nms = obj.groupNames{clrPos};
                    sz = length(nms);
                    for i = 1:sz
                        nm = nms {i};
                        cntr = cntr + 1;
                        plot(x, y, 'Color', lc{i}, 'LineStyle', '-', 'LineWidth', lthg);
                        hold on;
                        leg{cntr} = nm;
                    end        
                end

                if (lineStylePos > 0)
                    nms = obj.groupNames{lineStylePos};
                    sz = length(nms);
                    for i = 1:sz
                        nm = nms{i};
                        cntr = cntr + 1;
                        plot(x, y, 'Color', nolc, 'LineStyle', ls{i}, 'LineWidth', lthk(i));
                        hold on;
                        leg{cntr} = nm;
                    end
                end
                if (markerPos > 0)
                    nms = obj.groupNames{markerPos};
                    sz = length(nms);
                    for i = 1:sz
                        nm = nms {i};
                        cntr = cntr + 1;
                        plot(x, y, 'Color', clrwh, 'LineStyle', '-', 'LineWidth', 2, 'Marker', markers{i}, 'MarkerEdgeColor', nolc);
                        hold on;
                        leg{cntr} = nm;
                    end
                end
                if (option200_df_hack)
                        cntr = cntr + 1;
                        plot(x, y, 'Color', 'k', 'LineStyle', '-', 'LineWidth', 2);
                        hold on;
                        leg{cntr} = 'EPF';
                end
                if (option121_hack)
                   i = sz;
                   plot(x, y, 'Color', nolc, 'LineStyle', '-.', 'LineWidth', lthg);
                   cntr = cntr + 1;
                   leg{cntr} = '0D (EPF)';
                end
                lg = legend(leg, 'FontSize', obj.lfs, 'Interpreter','latex');
                legend('boxoff');


                % now plotting the data 
                for i = 1:numPts
                   ind = inds{i};
                   i1 = ind(1);
                   i2 = ind(2);
                   i3 = ind(3);
                   i4 = ind(4);
                   i5 = ind(5);
                   clr = [0, 0, 0];
                   lsv = '-';
                   lwv = lthg;
                   mrkr = 'none';
                   cntr = cntr + 1;
                   lc{cntr} = [0 0 0]; % black

                   if (~option200_df_hack || (i <= numPtsBeforeAdded))
                       if (clrPos > 0)
                           clr = lc{ind(clrPos)};
                       end
                       lstno = 1;
                       if (lineStylePos > 0)
                           lstno = ind(lineStylePos);
                       end
                       lsv = ls{lstno};
                       lwv = lthk(lstno);
                       if (markerPos > 0)
                           mrkr = markers{ind(markerPos)};
                       end
                   end
                   x = obj.xVals{i1}{i2}{i3}{i4}{i5};
                   dat = obj.data{i1}{i2}{i3}{i4}{i5};
                   y = dat(:, col);
                   if ((option121_hack) && (numPts - i < 3))
                       lsv = '-.';
                   end
                    plot(x, y, 'Color', clr, 'LineStyle', lsv, 'LineWidth', lwv, 'Marker', mrkr);
                    hold on;
                end
                fnbase = [folderName, '/', fnBase0, '_col_', icol];
                xh = get(gca, 'XLabel');
                set(xh, 'String', xLabel, 'FontSize', obj.labsz, 'VerticalAlignment','Top', 'Interpreter', 'latex');
                yh = get(gca, 'YLabel');
                set(yh, 'String', yLabel, 'FontSize', obj.labsz, 'VerticalAlignment','Bottom', 'Interpreter', 'latex');
            
                print('-dpng', [fnbase, '.png']);
                savefig([fnbase, '.fig']);
            end
            close('all');
        end
    end
end