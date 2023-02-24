classdef datawheader
    properties
        % angle-dependent computed values
        nFieldCols;
        fieldGroups;
        fieldSubGroups;
        fieldNames;
        fieldLatexNames;

        nDataRows;
        dataMat;
        rootFolder = '';
    end
    methods
        function objout = read(obj, fid)
            objout = obj;
            tline = fgetl(fid);
            C = strsplit(tline);
            objout.nFieldCols = length(C);
            for j = 1:objout.nFieldCols
                objout.fieldGroups{j} = C{j};
                objout.fieldSubGroups{j} = fscanf(fid, '%s', 1);
            end
            for j = 1:objout.nFieldCols
                objout.fieldNames{j} = fscanf(fid, '%s', 1);
            end
            for j = 1:objout.nFieldCols
                objout.fieldLatexNames{j} = fscanf(fid, '%s', 1);
            end
            objout.dataMat = fscanf(fid, '%g', [objout.nFieldCols, inf])';
            [objout.nDataRows, n] = size(objout.dataMat);
            for i = 1:objout.nDataRows
                for j = 1:n
                    tmp = objout.dataMat(i, j);
                    if (abs(tmp) > 0.99e40)
                        objout.dataMat(i, j) = nan;
                    end
                end
            end
        end
        function plotRaw(obj, startCol, prename)
            if nargin < 2
                startCol = 1;
            end
            if nargin < 3
                prename = '';
            end
            x = obj.dataMat(:,startCol);
            lfs = 18;
            for j = startCol + 1:obj.nFieldCols
                figure(j);
                y = obj.dataMat(:, j);
                plot(x, y);
                xh = get(gca, 'XLabel');
                set(xh, 'String', ['$$ ', obj.fieldLatexNames{startCol}, ' $$'], 'FontSize', lfs, 'VerticalAlignment','Top', 'Interpreter', 'latex');

                yh = get(gca, 'YLabel');
                set(yh, 'String', ['$$ ', obj.fieldLatexNames{j}, ' $$'], 'FontSize', lfs, 'VerticalAlignment','Bottom', 'Interpreter', 'latex');

                rt = '';
                if (strcmp(obj.rootFolder, '') == 0)
                    [status,message,messageid] = mkdir(obj.rootFolder);
                    rt = [obj.rootFolder, '/'];
                end
                name = [rt, prename, 'plot_', num2str(startCol, '%02d'), '_', num2str(j, '%02d'), '_', obj.fieldNames{j}];
                savefig([name, '.fig']);
                print('-dpng', [name, '.png']);
            end
        end
    end
end