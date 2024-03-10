function IO_shape_dd2(shape, dd2)
if (nargin < 1)
    shape = 1.5;
end
if (nargin < 2)
    dd2 = 0.5;
end
folder = ['dd', num2str(dd2)] ;
[status,msg,msgID] = mkdir(folder);
folder = [folder, '/shape', num2str(shape)];
[status,msg,msgID] = mkdir(folder);
folder = [folder, '/'];

if (nargin < 1)
    shape = 1.5;
end
if (nargin < 2)
    dd2 = 0.5;
end
llcs = -0.5:-0.5:-4.5;
nllcs = length(llcs);
for llci = 1:nllcs
    llc = llcs(llci);
    [sNs, mappedStats, snStats] = ReadStats(llc, shape, dd2);
    sz = length(sNs);
    for i = 1:sz
        sn = sNs(i);
        snStr = num2str(sn);
        mapstat = mappedStats(i, :);
        snstat = snStats(i, :);
        fn = [folder, snStr, '.txt'];
        fido = fopen(fn, 'a');
        fprintf(fido, '%g', llc);
        for j = 1:19
          fprintf(fido, '\t%g', mapstat(j));
        end
        for j = 1:19
          fprintf(fido, '\t%g', snstat(j));
        end
        fprintf(fido, '\n');
        fclose(fido);
    end
end
