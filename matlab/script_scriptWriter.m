versionEnd = 1619;
serialSt = 1000;
serialEn = 1039;
% isaac new gen: 54 for long_bigmem, 48 for campus
% isaac legacy: 24 for sigma
numProcessors = 24;

versionOffset = 0;
versionStart = 0;

for n = 1:numProcessors
    fn = ['script', num2str(n - 1), '.sh'];
    fidi{n} = fopen(fn, 'w');
end

cntr = 0;
for ser = serialSt:serialEn
    for vb = versionStart:versionEnd
        version = vb + versionOffset;
        res = mod(cntr, numProcessors);
        fid = fidi{res + 1};
        command = ['./crack -v ', num2str(version), ' -s ', num2str(ser)];
        fprintf(fid, command);
        fprintf(fid, '\n');
        cntr = cntr + 1;
    end
end

for n = 1:numProcessors
    fclose(fidi{n});
end
