function plot_InhomogeneousElasticFractureField_Results(fileBase)
if nargin < 1
    fileBase = 'TestInhomogeneousElasticFractorField';
end
fileName = [fileBase, '_index'];
dat1 = load([fileName, '.txt'], '-ascii');
[m, n] = size(dat1);
x = dat1(:, 1);
for i = 1:n - 1
    y = dat1(:, i + 1);
    figure(i);
    plot(x, y);
    print('-dpng', [fileName, num2str(i), '.png']);
end


fileName = [fileBase, '_x'];
dat2 = load([fileName, '.txt'], '-ascii');
[m, n] = size(dat2);
x = dat2(:, 1);
for i = 1:n - 1
    y = dat2(:, i + 1);
    figure(i + 10);
    plot(x, y);
    print('-dpng', [fileName, num2str(i), '.png']);
end
fclose('all');
close('all');