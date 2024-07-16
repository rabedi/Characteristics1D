clear all
close('all');
fclose('all');
%fid = fopen('config_periodic_color_plots.txt', 'r');
fid = fopen('config_periodic_color_plots_wLims.txt', 'r');
buf = fscanf(fid, '%s', 1);
ncol	= fscanf(fid, '%d', 1);
buf = fscanf(fid, '%s', 1);
pth	= fscanf(fid, '%s', 1);
buf = fscanf(fid, '%s', 1);
col_ys = fscanf(fid, '%d', 1);
buf = fscanf(fid, '%s', 1);
col_xlabel = fscanf(fid, '%s', 1);
buf = fscanf(fid, '%s', 1);
col_ylabel= fscanf(fid, '%s', 1);
nFiles = fscanf(fid, '%d', 1);
for fi = 1:nFiles
    buf = fscanf(fid, '%s', 1);
    fn{fi} = [pth, buf];
    xs(fi) = fscanf(fid, '%g', 1);
end
buf = fscanf(fid, '%s', 1);
sz = fscanf(fid, '%d', 1);
for i = 1:sz
    col(i) = fscanf(fid, '%d', 1);
    zmin = fscanf(fid, '%g', 1);
    zmax = fscanf(fid, '%g', 1);
    if (zmin > 1e39)
        zmin = nan;
    end
    if (zmax > 1e39)
        zmax = nan;
    end
    zmins(i) = zmin;
    zmaxs(i) = zmax;
end

% old data
%ncol = 108;
%ncol = 118;
%ncol = 126;
for fi = 1:nFiles
    fidr = fopen(fn{fi}, 'r');
    for j = 1:ncol
        header{j} = fscanf(fidr, '%s', 1);
    end
    cntr = 0;
    mat = [];
    while 1
        [vl, neof] = fscanf(fidr, '%f', 1);
        if (~neof)
            break;
        end
        cntr = cntr + 1;
        mat(cntr, 1) = vl;
        for j = 2:ncol
            mat(cntr, j) = fscanf(fidr, '%f', 1);
        end
    end
    fclose(fidr);
    ys{fi} = mat(:,col_ys);
    dat{fi} = mat;
end

xlbl = '$$\mathrm{log}_{10}(\bar{\dot{\epsilon}})$$';
ylbl = '$$\mathrm{log}_{10}(\bar{l})$$';

for i = 1:sz
    coln = col(i);
    zmin = zmins(i);
    zmax = zmaxs(i);
    colname = header{coln};
    for fi = 1:nFiles
        zs{fi} = dat{fi}(:,coln);
    end
    fnbase = ['Periodic_plot_', num2str(coln), colname];
    plot_xColumns_y_values_color(xs, ys, zs, xlbl, ylbl, fnbase, zmax, zmin);
    close('all');
    fclose('all');
end

