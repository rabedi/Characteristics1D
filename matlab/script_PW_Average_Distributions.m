clear all
delta = 0.5;
N = 100;
%alphas = [1.0, 1.5, 2.0, 3.0, 4.0];
alphas = [1.0, 1.5, 2.0, 4.0, 10.0, 20.0]; %, 100.0];
labsz = 22;
legsz = 14;

sz_alpha = length(alphas);
leg{1} = '$$ \alpha $$';
for ai = 1:sz_alpha
    alpha = alphas(ai);
    leg{ai + 1} = ['$$ ', num2str(alpha), ' $$'];
    [xs_tri{ai}, pdf_gen_tri{ai}, cdf_gen_tri{ai}, xs_Weibull{ai}, pdf_Weibull{ai}, cdf_Weibull{ai}] = GetPW_Average_Distributions(alpha, delta, N);
end
lc = cell(0);
cntr = 1;
lc{cntr} = [0 0 0]; % black
cntr = cntr + 1;
lc{cntr} = [1	0	0]; %red'
cntr = cntr + 1;
lc{cntr} = [1	102/255	0]; %orange'
cntr = cntr + 1;
%lc{cntr} = [0.5	0	1]; % purple 
%cntr = cntr + 1;
lc{cntr} = [0	0	1]; % blue 'b';
cntr = cntr + 1;
lc{cntr} = [0 135/255 0]; %green
cntr = cntr + 1;
lc{cntr} = [0	1	1]; % teal [0	0	0.5]; % dark blue
cntr = cntr + 1;
lc{cntr} = [0.5	0.25	0]; % brown
cntr = cntr + 1;
lc{cntr} = [0.5	0.5	0.5]; %'dark_gray2'
cntr = cntr + 1;
lc{cntr} = [0.75	0.75	0.75]; %'gray2'
cntr = cntr + 1;
lc{cntr} = [1	1	0]; %yellow [0.5	0.5	0]; % olive , [0	0.5	0.25]; %green blue
cntr = cntr + 1;

figure(1);
plot([nan], [nan], 'Color', 'w');
hold on;
for ai = 1:sz_alpha
    x = xs_tri{ai};
    y = pdf_gen_tri{ai};
    plot(x, y, 'LineWidth', 2, 'Color', lc{ai}); %, 'LineStyle', ls{j});
    hold on;
end
xh = get(gca, 'XLabel');
set(xh, 'String', '$$ s $$ (point-wise strength)', 'FontSize', labsz, 'VerticalAlignment','Top', 'Interpreter', 'latex');
yh = get(gca, 'YLabel');
set(yh, 'String', '$$ \mathrm{PDF} $$', 'FontSize', labsz, 'VerticalAlignment','Bottom', 'Interpreter', 'latex');
set(gca, 'ylim', [0, 15]);
legend(leg, 'FontSize', legsz, 'Interpreter', 'latex', 'Location', 'best');
legend('boxoff');
fn = 'zPDF_PW';
fne = [fn, '.png'];
print('-dpng', fne);
fne = [fn, '.fig'];
savefig(fne);


figure(2);
plot([nan], [nan], 'Color', 'w');
hold on;
for ai = 1:sz_alpha
    x = xs_tri{ai};
    y = cdf_gen_tri{ai};
    plot(x, y, 'LineWidth', 2, 'Color', lc{ai}); %, 'LineStyle', ls{j});
    hold on;
end
xh = get(gca, 'XLabel');
set(xh, 'String', '$$ s $$ (point-wise strength)', 'FontSize', labsz, 'VerticalAlignment','Top', 'Interpreter', 'latex');
yh = get(gca, 'YLabel');
set(yh, 'String', '$$ \mathrm{CDF} $$', 'FontSize', labsz, 'VerticalAlignment','Bottom', 'Interpreter', 'latex');
legend(leg, 'FontSize', legsz, 'Interpreter', 'latex', 'Location', 'best');
legend('boxoff');
fn = 'zCDF_PW';
fne = [fn, '.png'];
print('-dpng', fne);
fne = [fn, '.fig'];
savefig(fne);



%%%%%%%%%%%%%% Weibull's
figure(3);
plot([nan], [nan], 'Color', 'w');
hold on;
for ai = 1:sz_alpha
    x = xs_Weibull{ai};
    y = pdf_Weibull{ai};
    plot(x, y, 'LineWidth', 2, 'Color', lc{ai}); %, 'LineStyle', ls{j});
    hold on;
end
xh = get(gca, 'XLabel');
set(xh, 'String', '$$ S $$ (domain strength)', 'FontSize', labsz, 'VerticalAlignment','Top', 'Interpreter', 'latex');
yh = get(gca, 'YLabel');
set(yh, 'String', '$$ \mathrm{PDF} $$', 'FontSize', labsz, 'VerticalAlignment','Bottom', 'Interpreter', 'latex');
set(gca, 'ylim', [0, 30]);
legend(leg, 'FontSize', legsz, 'Interpreter', 'latex', 'Location', 'best');
legend('boxoff');
fn = 'zPDF_Domain';
fne = [fn, '.png'];
print('-dpng', fne);
fne = [fn, '.fig'];
savefig(fne);

figure(4);
plot([nan], [nan], 'Color', 'w');
hold on;
for ai = 1:sz_alpha
    x = xs_Weibull{ai};
    y = cdf_Weibull{ai};
    plot(x, y, 'LineWidth', 2, 'Color', lc{ai}); %, 'LineStyle', ls{j});
    hold on;
end
xh = get(gca, 'XLabel');
set(xh, 'String', '$$ S $$ (domain strength)', 'FontSize', labsz, 'VerticalAlignment','Top', 'Interpreter', 'latex');
yh = get(gca, 'YLabel');
set(yh, 'String', '$$ \mathrm{CDF} $$', 'FontSize', labsz, 'VerticalAlignment','Bottom', 'Interpreter', 'latex');
legend(leg, 'FontSize', legsz, 'Interpreter', 'latex', 'Location', 'best');
legend('boxoff');
fn = 'zCDF_Domain';
fne = [fn, '.png'];
print('-dpng', fne);
fne = [fn, '.fig'];
savefig(fne);



%%%%%%%%%%
Ns = [10, 100, 1000, 10000];
sz_N = length(Ns);
legn{1} = '$$ N $$';
alpha = 2.0;
for ni = 1:sz_N
    N = Ns(ni);
    legn{ni + 1} = ['$$ ', num2str(N), ' $$'];
    [xs_trin{ni}, pdf_gen_trin{ni}, cdf_gen_trin{ni}, xs_Weibulln{ni}, pdf_Weibulln{ni}, cdf_Weibulln{ni}] = GetPW_Average_Distributions(alpha, delta, N);
end


%%%%%%%%%%%%%% Weibull's - effect of N
figure(5);
plot([nan], [nan], 'Color', 'w');
hold on;
for ni = 1:sz_N
    x = xs_Weibulln{ni};
    y = pdf_Weibulln{ni};
    plot(x, y, 'LineWidth', 2, 'Color', lc{ni + 5}); %, 'LineStyle', ls{j});
    hold on;
end
xh = get(gca, 'XLabel');
set(xh, 'String', '$$ S $$ (domain strength)', 'FontSize', labsz, 'VerticalAlignment','Top', 'Interpreter', 'latex');
yh = get(gca, 'YLabel');
set(yh, 'String', '$$ \mathrm{PDF} $$', 'FontSize', labsz, 'VerticalAlignment','Bottom', 'Interpreter', 'latex');
set(gca, 'ylim', [0, 30]);
legend(legn, 'FontSize', legsz, 'Interpreter', 'latex', 'Location', 'best');
legend('boxoff');
fn = 'zPDF_Domain_NEffect';
fne = [fn, '.png'];
print('-dpng', fne);
fne = [fn, '.fig'];
savefig(fne);

figure(6);
plot([nan], [nan], 'Color', 'w');
hold on;
for ni = 1:sz_N
    x = xs_Weibulln{ni};
    y = cdf_Weibulln{ni};
    plot(x, y, 'LineWidth', 2, 'Color', lc{ni + 5}); %, 'LineStyle', ls{j});
    hold on;
end
xh = get(gca, 'XLabel');
set(xh, 'String', '$$ S $$ (domain strength)', 'FontSize', labsz, 'VerticalAlignment','Top', 'Interpreter', 'latex');
yh = get(gca, 'YLabel');
set(yh, 'String', '$$ \mathrm{CDF} $$', 'FontSize', labsz, 'VerticalAlignment','Bottom', 'Interpreter', 'latex');
legend(legn, 'FontSize', legsz, 'Interpreter', 'latex', 'Location', 'best');
legend('boxoff');
fn = 'zCDF_Domain_NEffect';
fne = [fn, '.png'];
print('-dpng', fne);
fne = [fn, '.fig'];
savefig(fne);
