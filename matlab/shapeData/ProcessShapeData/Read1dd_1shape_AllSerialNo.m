function [dat_llc_indexed, serNo_llc_indexed, cntrs_llc_indexed] = Read1dd_1shape_AllSerialNo(shapeStr, ddStr, serNoMax)
if (nargin < 1)
    shapeStr = '2';
end
if (nargin < 2)
    ddStr = '0.9';
end
if (nargin < 3)
    serNoMax = 4999;
end
fprintf(1, '\tdd = %s, shape = %s\n', ddStr, shapeStr);
dat_llc_indexed = cell(9, 1);
serNo_llc_indexed = cell(9, 1);
cntrs_llc_indexed = zeros(9, 1);

for serNo = 1:serNoMax
    [dataOut, exists, validFile] = Read1dd_1shape_1ser(serNo, shapeStr, ddStr);
    if (validFile == 0)
        break;
    end
    for j = 1:9
        if (exists(j) == 1)
            cntrs_llc_indexed(j) = cntrs_llc_indexed(j) + 1;
            cntr = cntrs_llc_indexed(j);
            serNo_llc_indexed{j}(cntr) = serNo;
            dat_llc_indexed{j}(cntr,:) = dataOut{j};
            n = length(dataOut{j});
        end
    end
end
[status,msg,msgID] = mkdir('stat');
[status,msg,msgID] = mkdir('pdf');
[status,msg,msgID] = mkdir('Weibull');

for fld = 1:n
    [status,msg,msgID] = mkdir(['stat/fld', num2str(fld)]);
    [status,msg,msgID] = mkdir(['pdf/fld', num2str(fld)]);
    [status,msg,msgID] = mkdir(['Weibull/fld', num2str(fld)]);
end

file_out_base = ['shape', shapeStr, '_dd', ddStr];
for j = 1:9
    corLen = -0.5 + (j - 1) * -0.5;
    corLenStr = num2str(corLen);
    fileBase = [file_out_base, '_llc', corLenStr];
    for fld = 1:n
        fldStr = num2str(fld);
        fileNameBase2 = ['fld', fldStr, '/', fileBase, '_fld', num2str(fld)];
        vals = dat_llc_indexed{j}(:,fld);
        opss = OnePointStatSimple;
        min4Weibull = 1 - str2num(ddStr);
        doWeibull = (fld == 2);
        doGauss = 1;
        opss = opss.ComputeVals(vals, min4Weibull, doWeibull, doGauss);
        % plot(opss.pdf_xAxis, opss.pdf_yAxis);
        fid = fopen(['stat/', fileNameBase2, '_stat.txt'], 'w');
        printHeader = 1;
        opss.printVals(fid, printHeader);
        fclose(fid);
        fid = fopen(['pdf/', fileNameBase2, '_pdf.txt'], 'w');
        sz = length(opss.pdf_xAxis);
        for k = 1:sz
            fprintf(fid, '%g\t', opss.pdf_xAxis(k));
        end
        fprintf(fid, '\n');
        for k = 1:sz
            fprintf(fid, '%g\t', opss.pdf_yAxis(k));
        end
        fprintf(fid, '\n');
        fclose(fid);

        if (doWeibull)
            % Create the structure
            WeibullR_struct = struct('p', opss.WeibullR_p, ...
                         'xreg', opss.WeibullR_xreg, ...
                         'yreg', opss.WeibullR_yreg, ...
                         'X', opss.WeibullR_X, ...
                         'Y', opss.WeibullR_Y, ...
                         'R2', opss.WeibullR_R2, ...
                         'shape', opss.WeibullR_shape, ...
                         'scale', opss.WeibullR_scale, ...
                         'N', opss.WeibullR_N);

            % Save the structure to a binary file
            fileName = ['Weibull/', fileNameBase2, '_Weibull.mat'];
            save(fileName, 'WeibullR_struct', '-v7.3'); % Use -v7.3 for large data
        end
    end
end
