function plot_InhomogeneousField_Results(fileBase)
if nargin < 1
    fileBase = 'file_v1_sz0_1';
    fileBase = 'file_v0_sz0_1';
%    fileBase = 'file_v1_sz0_x_1';
    fileBase = 'file_v0_sz0_x_1';
%    fileBase = 'file0';
%    fileBase = 'file_nox_0';
end
dat1 = load([fileBase, '_vals.txt'], '-ascii');
figure(1);
plot(dat1);
print('-dpng', [fileBase, '_vals.png']);

dat2 = load([fileBase, '_x_v.txt'], '-ascii');
x = dat2(:, 1);
y = dat2(:, 2);
figure(2);
plot(x, y);
print('-dpng', [fileBase, '_x_v.png']);
