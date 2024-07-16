function plot_xColumns_y_values_color(xs, ys, zs, xlbl, ylbl, fnbase, zmax, zmin)
labsz = 25;

if (nargin < 1)
    xs = [0, 1, 2];
end
if (nargin < 2)
    ys{1} = [0, 1, 2, 3, 4, 5];
    ys{2} = [0.5, 1.5, 2.5, 3.5];
    ys{3} = [0, 1, 2, 3, 3.5, 4];
end
if (nargin < 3)
    zs{1} = [13, 2, 1, 0, 1, 2];
    zs{2} = [1, 1, 1, 1];
    zs{3} = [0, 1, 2, 3, 4, 15];
end    
if (nargin < 4)
    xlbl = '$$ \bar{\dot{\epsilon}} $$';
end
if (nargin < 5)
    ylbl = 'ylabel';
end
if (nargin < 6)
    fnbase = 'colorplot_example';
end
if (nargin < 7)
    zmax = nan;
    zmax = 5;
end
if (nargin < 8)
    zmin = nan;
    zmin = 0;
end

    

% xs: vector of xs, e.g. xs(1), ...
% nx = length(xs)
% ys: cell of size nx, each ys{i} contains some y values
% zs: same as ys, but contains z values

nx = length(xs);
nys = zeros(nx, 1);
for i = 1:nx
    nys(i) = length(ys{i});
end
if (isnan(zmax))
    zmax = -inf;
    zmin = inf;
    for i = 1:nx
        zmax = max(zmax, max(zs{i}));
        zmin = min(zmin, min(zs{i}));
    end
end
zFactor = 1.0 / (zmax - zmin);

% Choose a colormap (e.g., 'jet', 'parula', etc.)
cmap = jet; % You can also try 'parula', 'hot', 'cool', etc.
%cmap = parula; % You can also try 'parula', 'hot', 'cool', etc.


% Map the normalized z values to colors using the chosen colormap
num_colors = size(cmap, 1);

% Create a figure
figure;

for i = 1:nx - 1
    ip = i;
    in = i + 1;
    nyp = nys(ip);
    nyn = nys(in);
    nmin = min(nyp, nyn);
    cntr = 0;
    for j = 1:nmin - 1
        for k = 1:2
            if (k == 1)
                xindices = [ip, in, ip];
                yindices = [j, j, j + 1];
            else
                xindices = [ip, in, in];
                yindices = [j + 1, j, j + 1];
            end
            for m = 1:3
                xind = xindices(m);
                vertices(m, 1) = xs(xind);
                vertices(m, 2) = ys{xind}(yindices(m));
                z_values(m) = zs{xind}(yindices(m));
                z_normalized(m) = zFactor * (z_values(m) - zmin);
                if (z_normalized(m) > 1.0)
                    z_normalized(m) = 1.0;
                elseif (z_normalized(m) < 0.0)
                    z_normalized(m) = 0.0;
                end
                % z_normalized(m) = zmin + (zmax - zmin) * z_normalized(m);
            end
            colors = interp1(linspace(0, 1, num_colors), cmap, z_normalized);
            % Plot the triangle with interpolated colors
            patch('Vertices', vertices, 'Faces', [1, 2, 3], ...
                  'FaceVertexCData', colors, ...
                  'FaceColor', 'interp', ...
                  'EdgeColor', 'none'); % Remove edge color for better interpolation
        end
        if (nyp == nyn)
            continue;
        end
        xindices_set = cell(0);
        yindices_set = cell(0);
        cntr_t = 0;
        if (nyp < nyn)
            for k = nyp:nyn - 1
                xindices = [ip, in, in];
                yindices = [nyp, k, k + 1];
                cntr_t = cntr_t + 1;
                xindices_set{cntr_t} = xindices;
                yindices_set{cntr_t} = yindices;
            end
        else
            for k = nyn:nyp - 1
                xindices = [ip, in, ip];
                yindices = [k, nyn, k + 1];
                cntr_t = cntr_t + 1;
                xindices_set{cntr_t} = xindices;
                yindices_set{cntr_t} = yindices;
            end
        end
        for k = 1:cntr_t
            xindices = xindices_set{k};
            yindices = yindices_set{k};

            for m = 1:3
                xind = xindices(m);
                vertices(m, 1) = xs(xind);
                vertices(m, 2) = ys{xind}(yindices(m));
                z_values(m) = zs{xind}(yindices(m));
                z_normalized(m) = zFactor * (z_values(m) - zmin);
                if (z_normalized(m) > 1.0)
                    z_normalized(m) = 1.0;
                elseif (z_normalized(m) < 0.0)
                    z_normalized(m) = 0.0;
                end
                % z_normalized(m) = zmin + (zmax - zmin) * z_normalized(m);
            end
            colors = interp1(linspace(0, 1, num_colors), cmap, z_normalized);
            % Plot the triangle with interpolated colors
            patch('Vertices', vertices, 'Faces', [1, 2, 3], ...
                  'FaceVertexCData', colors, ...
                  'FaceColor', 'interp', ...
                  'EdgeColor', 'none'); % Remove edge color for better interpolation
        end
    end
end

% Set the axis limits and equal aspect ratio
%axis([0 1 0 1]);
%axis equal;
%grid on;

%title('Triangle with Interpolated Colors Based on Z Values');
cmap = jet; % You can also try 'parula', 'hot', 'cool', etc.
colormap(cmap);
caxis([zmin, zmax]);
colorbar; % Add a colorbar to show the mapping of z values to colors

xh = get(gca, 'XLabel');
set(xh, 'String', xlbl, 'FontSize', labsz, 'VerticalAlignment','Top', 'Interpreter', 'latex');
yh = get(gca, 'YLabel');
set(yh, 'String', ylbl, 'FontSize', labsz, 'VerticalAlignment','Bottom', 'Interpreter', 'latex');

print('-dpng', [fnbase, '.png']);
savefig([fnbase, '.fig']);
