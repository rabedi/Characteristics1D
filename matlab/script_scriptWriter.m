versionEnd = 1619;
serialSt = 1200;
serialEn = 1239;
% isaac new gen: 54 for long_bigmem, 48 for campus
% isaac legacy: 24 for sigma
numProcessors = 48;
versionOffset = 0;
versionStart = 0;


% la -2
versionStart = 1620;
versionEnd = 1754;
% la -2p5
versionStart = 1755;
versionEnd = 1889;
% la -3
versionStart = 1890;
versionEnd = 2024;
numProcessors = 56;

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
