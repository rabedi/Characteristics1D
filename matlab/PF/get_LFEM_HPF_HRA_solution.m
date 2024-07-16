function [eps, D, sig, sigmap_Max_eps, sigmap_Max, eps_f, phi, phi_unloading, phi_loading, brittleness_phi, brittleness_strain] = get_LFEM_HPF_HRA_solution(model_s, CZM_normalization4AT1_2, doPlot)

if (nargin < 1)
    model_s = 'AT1';
%    model_s = 'AT2';
end
if (nargin < 2)
    % 0: no normalization: lScale = b
    % 1: quasi-static strength = 1 
    % 2: A is made equal to 1: useful for HRA-HPF
    CZM_normalization4AT1_2 = 0;
    CZM_normalization4AT1_2 = 1;
%    CZM_normalization4AT1_2 = 2;
end

if (nargin < 3)
    doPlot = 1;
end

[F, A, beta] = get_LFEM_PH_F_A_beta(model_s, CZM_normalization4AT1_2);
del_eps = 0.00001;
eps_f = 6;
eps = 0:del_eps:eps_f;
eps(1) = 0.000001 * eps_f;
nu = 0.25;

sqrt_x = sqrt(eps);
Ax2 = A .* eps .* eps;
c2 = gamma(1 - nu) * power(0.5 * A, nu);
omD = c2 .* sqrt_x .* besselj(-nu, Ax2);
ind = find(omD < 0);
sz = ind(1) - 1;
eps = eps(1:sz);
omD = omD(1:sz);
D = 1.0 - omD;
sigma = omD .* omD .* eps;

[sigmap_Max, i] = max(sigma);
sigmap_Max_eps = eps(i);
eps_f = eps(sz);
phi = trapz(eps, sigma);
phi_unloading = trapz(eps(i:sz), sigma(i:sz));
phi_loading = phi - phi_unloading;
brittleness_phi = phi_loading / phi;
brittleness_strain = sigmap_Max_eps / eps_f;


if (doPlot)
    format long;
    filenameBase = ['HPF_HRA_', model_s, '_CZM_normalization_', num2str(CZM_normalization4AT1_2)];
    fid = fopen([filenameBase, '.txt'], 'w');
    fprintf(fid, 'epsM %f\n', sigmap_Max_eps);
    fprintf(fid, 'sigM %f\n', sigmap_Max);
    fprintf(fid, 'phiM %f\n', phi_loading);
    fprintf(fid, 'epsF %f\n', eps_f);
    fprintf(fid, 'phiF %f\n', phi);
    fprintf(fid, 'brittleness_phi %f\n', brittleness_phi);
    fprintf(fid, 'phi_unloading %f\n', phi_unloading);
    fprintf(fid, 'brittleness_strain %f\n', brittleness_strain);
    fclose(fid);

    clf;
    sigmap_Max
    sigmap_Max_eps 
    eps_f
    phi
    phi_unloading
    phi_loading
    brittleness_phi
    brittleness_strain

%    figure(1);
%    plot(eps, omD);
    
    figure(1);
    plot(eps, sigma, '-r');
    hold on;
    plot(eps, D, '-b');
    lg = legend({'$$ \sigma $$', '$$ D$$ '}, 'FontSize', 18, 'Interpreter', 'latex', 'Location', 'NorthWest');
    legend('boxoff');

    xlab = '$$ \hat{\bar{\epsilon}} $$';
    ylab = 'Ylabel';
    labsz = 24;
    xh = get(gca, 'XLabel');
    set(xh, 'String', xlab, 'FontSize', labsz, 'VerticalAlignment','Top', 'Interpreter', 'latex');
    yh = get(gca, 'YLabel');
    set(yh, 'String', ylab, 'FontSize', labsz, 'VerticalAlignment','Bottom', 'Interpreter', 'latex');

    print('-dpng', [filenameBase, '.png']);
    savefig([filenameBase, '.fig']);
end