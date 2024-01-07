function [epsilon_p, sigma_p, D_vec, Dp_vec, Dpp_vec, omegaD_vec, sigmap_Max, sigmap_Max_eps, eps_f, ...
    phi, phi_unloading, phi_loading, brittleness_phi, brittleness_strain] = ...
    ComputePF_withLoadingRate(ap, useAnalyticalSln, isHyper, xi, omegaCZM, bPrime, CZM_model_name, CZM_normalization4AT1_2)
%CZM_normalization4AT1_2: means that when we use AT1/AT2 models, we change
%the normalization from sigmaNormalization = sqrt(G/Eb) and
%lengthNormalization = b to sigmaNormalization = sigma_coh, and
%lengthNormalization = l_coh. This option = 1 is useful if we compare AT1
%and AT2 with CZM models

%eps, sigma, tpp, 
delD = 1e-3;
dely = 1e-4;
% Newmark
delta = 0.5;

if nargin < 1
    ap = 100;
end

if nargin < 2
    useAnalyticalSln = 0;
end
if nargin < 3
    isHyper = 1;
end
if nargin < 4
    xi = 1;
end
if nargin < 5
    if (xi == 2)
        omegaCZM = 1;
    else
        omegaCZM = 0;
    end
end

if nargin < 6
    bPrime = -1.0;
end

if nargin < 7
    CZM_model_name = 'Linear';
end

if nargin < 8
    CZM_normalization4AT1_2 = 1;
end

[yEpsp_vec, sigmap_vec, D_vec, Dp_vec,  Dpp_vec, omegaD_vec, sigmap_Max, sigmap_Max_y, yEpsp_f, ...
    phi, phi_unloading, phi_loading, brittleness_phi, brittleness_strain] = ...
    ComputePF(useAnalyticalSln, isHyper, xi, omegaCZM, bPrime, CZM_model_name, CZM_normalization4AT1_2);

% epsilon_p = epsilon_ini_p + y * factor
% sigma_p = omega(D) * epsilon_p
% epsilon = 
% ap^(-1/2) H, ap^(-2/3) P
% factor = 1/epsilon for H, 1/sqrt(epsilon) for P
% = sqrt(a) H, a^(1/3)
if (isHyper)
    factor = sqrt(ap);
else
    factor = power(ap, 1.0/3.0);
end
epsilon_p_ini = 1.0; % initial strength - correct for CZM
if (xi == 0) % AT2
    epsilon_p_ini = 0.0;
elseif (xi == 1) % AT1
    if (CZM_normalization4AT1_2 == 0) % CZM normalization is not used, so the limit of linear response is sqrt(3/8)
        epsilon_p_ini = sqrt(3.0/8.0);
    end % otherwise is 1.0
end
epsilon_p = epsilon_p_ini + factor * yEpsp_vec;
sigma_p = omegaD_vec .* epsilon_p;

phi_ini = 0.5 * epsilon_p_ini * epsilon_p_ini;
sz = length(epsilon_p);
[sigmap_Max, i] = max(sigma_p);
sigmap_Max_eps = epsilon_p(i);
eps_f = epsilon_p(sz);
epsilon_p_diff = epsilon_p - epsilon_p_ini;
phi = trapz(epsilon_p_diff, sigma_p) + phi_ini;
phi_unloading = trapz(epsilon_p_diff(i:sz), sigma_p(i:sz));
phi_loading = phi - phi_unloading;
brittleness_phi = phi_loading / phi;
brittleness_strain = sigmap_Max_eps / eps_f;

epsilon_p = [[0] epsilon_p];
sigma_p = [[0] sigma_p];
D_vec = [[0] D_vec];
omegaD_vec = [[1] omegaD_vec];
Dp_vec = [[0] Dp_vec];
if (length(Dpp_vec) > 0)
    Dpp_vec = [[0] Dpp_vec];
end

if (1)
    figure(1);
    plot(epsilon_p, sigma_p);
    figure(2);
    plot(epsilon_p, D_vec);
    figure(3);
    plot(epsilon_p, omegaD_vec);
    a = 12;
end
