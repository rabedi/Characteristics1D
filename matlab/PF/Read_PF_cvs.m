function [data, labels] = Read_PF_cvs(fileName)
if (nargin < 1)
    fileName = 'option_31_H_AT1_lcRat_none_df_none/scalars.csv';
end

if (contains(fileName, 'model') == 0)
    %M = csvread(fileName);
    M = readmatrix(fileName);
    [m, n] = size(M);
    data = M(:, 11:n);
    labels = M(:,1)';
else
    fileName_new = 'option_5_ap_-1_isHyper_-1/scalars.csv';
    M = readmatrix(fileName_new);
    [m, n] = size(M);
    row = 0;
    if (contains(fileName, 'AT1'))
        row = 1;
    elseif (contains(fileName, 'AT2'))
        row = 2;
    elseif (contains(fileName, 'CZM-W'))
        row = 3;
    elseif (contains(fileName, 'CZM-L'))
        row = 4;
    end
    datShort = M(row, 2:9);
    labels = -4:0.25:6;
    m = length(labels);
    data = zeros(m, 8);
    for i = 1:m
        for j = 1:8
            data(i, j) = datShort(j);
        end
    end
end

jis = 1:6;
for i = 1:m
    for ji = 1:6
        j = jis(ji);
        data(i, j) = log10(data(i, j));
    end
end

