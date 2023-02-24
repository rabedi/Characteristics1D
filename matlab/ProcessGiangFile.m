function ProcessGiangFile(serNum)
if nargin < 1
    serNum = 984;
end
for cl = -1:-0.5:-1
    clstr = num2str(cl, '%0.1f');
    for serNum = 1:5000
        sers = num2str(serNum);
        fileName = ['../InhomogeneousFiles/data_n10_cl', clstr, '_np1025/initial_values', sers, '_n10.txt'];
        sersB = num2str(serNum - 1);
        fnout = ['../InhomogeneousFiles/cl', clstr, '_np1025/initial_values_', sersB, '.txt'];
    %    copyfile(fileName, fnout);
    
        fid = fopen(fileName, 'r');
        if (fid < 0)
            fileName
        end
        
        num = fscanf(fid, '%d', 1);
        dat = fscanf(fid, '%f', inf);
        fclose(fid);
        %[f,xi] = ksdensity(dat);
        %meanV = mean(dat);
        %sdiv = std(dat);
        %plot(xi,f);
        
        if (serNum == -1)
            x = 0:1/(num - 1):1;
            plot(x, dat)
        end
        fid = fopen(fnout, 'w');
        for i = 1:num
            fprintf(fid, '%g\t', dat(i));
        end
        fclose(fid);
        fprintf(1, '%d\t', serNum);
    end
    fprintf(1, '\n%g finished\n', cl);
end