isOrtiz = 0;
plot_laMode = 0; % 0 actual a, 1 a*2 is used, 2 a/G is used 
add_dilute_approx = 0; %0;
add_Grady = 1;
add_Glenn = 1;
add_max_ene_real = 1;
add_max_ene_max = 0; %0;
add_vs0Final = 1; %0; v0 is zero at failure stage
add_delus0 = 1; %0;
add_delusm1 = 0; %0;
add_vs0Max = 1; %0; % v0 is zero at max stage -> only for intrinsic models -> Drugan model
if (isOrtiz == 1)
    add_vs0Max = 0;
end

labsz = 22;
legsz = 14;

cntr = 0;
leg = cell(0);
ls = cell(0);
lc = cell(0);
dat = cell(0);
num_pts = cell(0);
la = cell(0);
la2p = cell(0);
laGp = cell(0);

cntr = cntr + 1;
leg{cntr} = '\mathrm{max}(\mathcal{E}(t_{\mathrm{dilute}}))';
ls{cntr} = '-';
lc{cntr} = [0 0 0]; % black
[dat{cntr}, names, num_pts{cntr}, la{cntr}, la2p{cntr}, laGp{cntr}] = Periodic_Segment_ReadOneSolve('../_Solver1D_pfragment__SlnCntr_2_AddedPts0_solveNum_0.txt');


cntr = cntr + 1;
leg{cntr} = '\mathrm{dilute}';
ls{cntr} = '-';
lc{cntr} = [0.5	0.5	0.5]; %'dark_gray2'
[dat{cntr}, names, num_pts{cntr}, la{cntr}, la2p{cntr}, laGp{cntr}] = Periodic_Segment_ReadOneSolve('../_Solver1D_pfragment_PrimaryPt_0.txt');

if add_dilute_approx
    cntr = cntr + 1;
    leg{cntr} = '\mathrm{dilute}_\mathrm{approx}';
    ls{cntr} = '-';
    lc{cntr} = [0.75	0.75	0.75]; %'gray2'
    [dat{cntr}, names, num_pts{cntr}, la{cntr}, la2p{cntr}, laGp{cntr}] = Periodic_Segment_ReadOneSolve('../_Solver1D_pfragment_PrimaryPt_1.txt');
end

cntr = cntr + 1;
leg{cntr} = '\mathrm{Zhu-a}';
ls{cntr} = '--';
lc{cntr} = [0	0	1]; % blue 'b';
[dat{cntr}, names, num_pts{cntr}, la{cntr}, la2p{cntr}, laGp{cntr}] = Periodic_Segment_ReadOneSolve('../_Solver1D_pfragment_PrimaryPt_2.txt');

cntr = cntr + 1;
leg{cntr} = '\mathrm{Zhu-b}';
ls{cntr} = '--';
lc{cntr} = [1	0	0]; %red'
[dat{cntr}, names, num_pts{cntr}, la{cntr}, la2p{cntr}, laGp{cntr}] = Periodic_Segment_ReadOneSolve('../_Solver1D_pfragment_PrimaryPt_3.txt');


if add_Glenn
    cntr = cntr + 1;
    leg{cntr} = '\mathrm{Glenn}';
    ls{cntr} = '-.';
    lc{cntr} = [0	1	1]; % teal [0	0	0.5]; % dark blue
    [dat{cntr}, names, num_pts{cntr}, la{cntr}, la2p{cntr}, laGp{cntr}] = Periodic_Segment_ReadOneSolve('../_Solver1D_pfragment_PrimaryPt_4.txt');
end

if add_Grady
    cntr = cntr + 1;
    leg{cntr} = '\mathrm{Grady}';
    ls{cntr} = '-.';
    lc{cntr} = [1	102/255	0]; %orange'
    [dat{cntr}, names, num_pts{cntr}, la{cntr}, la2p{cntr}, laGp{cntr}] = Periodic_Segment_ReadOneSolve('../_Solver1D_pfragment_PrimaryPt_5.txt');
end

if (add_max_ene_real)
    cntr = cntr + 1;
    ls{cntr} = '-';
    lc{cntr} = [0 135/255 0]; %green
    leg{cntr} = '\mathrm{max}(\mathcal{E}/{\mathcal{E}_{\mathrm{src}}})';
    [dat{cntr}, names, num_pts{cntr}, la{cntr}, la2p{cntr}, laGp{cntr}] = Periodic_Segment_ReadOneSolve('../_Solver1D_pfragment__SlnCntr_2_AddedPts0_solveNum_1.txt');
end

if (add_max_ene_max)
    cntr = cntr + 1;
    ls{cntr} = '-';
    lc{cntr} = [1	1	0]; %yellow [0.5	0.5	0]; % olive , [0	0.5	0.25]; %green blue
    leg{cntr} = '\mathrm{max}(\mathcal{E}/{\mathcal{E}_{\mathrm{src,max}}})';
    [dat{cntr}, names, num_pts{cntr}, la{cntr}, la2p{cntr}, laGp{cntr}] = Periodic_Segment_ReadOneSolve('../_Solver1D_pfragment__SlnCntr_2_AddedPts0_solveNum_2.txt');
end

if (add_vs0Max)
    cntr = cntr + 1;
    ls{cntr} = ':';
    lc{cntr} = [0.5	0	1]; % purple 
    leg{cntr} = 'v_s^m = 0';
    [dat{cntr}, names, num_pts{cntr}, la{cntr}, la2p{cntr}, laGp{cntr}] = Periodic_Segment_ReadOneSolve('../_Solver1D_pfragment__SlnCntr_2_AddedPts0_solveNum_6.txt');
end

if (add_vs0Final)
    cntr = cntr + 1;
    ls{cntr} = ':';
    lc{cntr} = [0.5	0.25	0]; % brown
    leg{cntr} = 'v_s^f = 0';
    [dat{cntr}, names, num_pts{cntr}, la{cntr}, la2p{cntr}, laGp{cntr}] = Periodic_Segment_ReadOneSolve('../_Solver1D_pfragment__SlnCntr_2_AddedPts0_solveNum_3.txt');
end

if (add_delus0)
    cntr = cntr + 1;
    ls{cntr} = ':';
    lc{cntr} = [0.5	0	0.25]; %arghavani [1	0.5	0.5]; % peach
    leg{cntr} = '{\bar{\epsilon}}_s = 0';
    [dat{cntr}, names, num_pts{cntr}, la{cntr}, la2p{cntr}, laGp{cntr}] = Periodic_Segment_ReadOneSolve('../_Solver1D_pfragment__SlnCntr_2_AddedPts0_solveNum_4.txt');
end

if (add_delusm1)
    cntr = cntr + 1;
    ls{cntr} = ':';
    lc{cntr} = [1	0	1]; % magenta
    leg{cntr} = '{\bar{\epsilon}}_s = -1';
    [dat{cntr}, names, num_pts{cntr}, la{cntr}, la2p{cntr}, laGp{cntr}] = Periodic_Segment_ReadOneSolve('../_Solver1D_pfragment__SlnCntr_2_AddedPts0_solveNum_5.txt');
end

sz = cntr;
xlab = '$$ log_{10}(\dot{\epsilon}/\tilde{\dot{\epsilon}}) $$';
ylab = '$$ y $$';
for j = 1:sz
    leg{j} = ['$$ ', leg{j}, ' $$'];
end

num_y = length(names);
for yi = 1:num_y
    figure(1);
    ys = num2str(yi);
    fn = ['y-', ys, '-', names{yi}];
    
    for j = 1:sz
        x = la{j};
        if (plot_laMode == 1)
            x = la2p{j};
        elseif (plot_laMode == 2)
            x = laGp{j};
        end
        y = dat{j}(:, yi);
        plot(x, y, 'LineWidth', 2, 'Color', lc{j}, 'LineStyle', ls{j});
        hold on;
    end
    xh = get(gca, 'XLabel');
    set(xh, 'String', xlab, 'FontSize', labsz, 'VerticalAlignment','Top', 'Interpreter', 'latex');
    yh = get(gca, 'YLabel');
    set(yh, 'String', ylab, 'FontSize', labsz, 'VerticalAlignment','Bottom', 'Interpreter', 'latex');
    legend(leg, 'FontSize', legsz, 'Interpreter', 'latex', 'Location', 'best');
    legend('boxoff');

    fne = [fn, '.png'];
    print('-dpng', fne);
    fne = [fn, '.fig'];
    savefig(fne);
    clf;
end
close('all')
fclose('all')