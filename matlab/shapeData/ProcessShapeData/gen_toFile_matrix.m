function gen_toFile_matrix(fido, mat)
[m, n] = size(mat);
if ((m == 1) && (n == 1))
    fprintf(fido, '%g\t', mat(1, 1));
    return;
end
if ((m == 0) && (n == 0))
    fprintf(fido, 'null\t');
    return;
end
prntRow = ((m > 3) && (n > 3));
fprintf(fido, '[ ( %d , %d )', m, n);
for i = 1:m
    if (prntRow)
        fprintf(fido, '\t[row%d]', i);
    end
    for j = 1:n
        fprintf(fido, '\t%g', mat(i, j));
    end
end
fprintf(fido, ' ]\t\n');



