function mat = gen_fromFile_matrix(fidi)
buf = READ(fidi, 's');
mat = str2num(buf);
if (length(mat) == 1)
    return;
end
if (strcmp(buf, 'null') == 1)
    mat = [];
    return;
end

buf = READ(fidi, 's');
m = READ(fidi, 'd');
buf = READ(fidi, 's');
n = READ(fidi, 'd');
mat = zeros(m, n);

prntRow = ((m > 3) && (n > 3));

buf = fscanf(fidi, '%s', 1);
for i = 1:m
    if (prntRow)
        buf = fscanf(fidi, '%s', 1);
    end
    for j = 1:n
        mat(i, j) = fscanf(fidi, '%g', 1);
    end
end
buf = fscanf(fidi, '%s', 1);
