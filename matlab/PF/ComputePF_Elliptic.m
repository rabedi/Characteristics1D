function [epsilon_p_vec, sigma_p_vec, D_vec, Dp_vec, Dpp_vec, omegaD_vec, sigmap_Max, sigmap_Max_eps, eps_f, ...
    phi, phi_unloading, phi_loading, brittleness_phi, brittleness_strain] = ...
    ComputePF_Elliptic(xi, omegaCZM, bPrime, CZM_model_name, CZM_normalization4AT1_2)
%CZM_normalization4AT1_2: means that when we use AT1/AT2 models, we change
%the normalization from sigmaNormalization = sqrt(G/Eb) and
%lengthNormalization = b to sigmaNormalization = sigma_coh, and
%lengthNormalization = l_coh. This option = 1 is useful if we compare AT1
%and AT2 with CZM models

%eps, sigma, tpp, 
delD = 1e-3;

if nargin < 1
    xi = 0;
end
if nargin < 2
    if (xi == 2)
        omegaCZM = 1;
    else
        omegaCZM = 0;
    end
end

if nargin < 3
    bPrime = -1;
end

if nargin < 4
    CZM_model_name = 'Linear';
end

if nargin < 5
    CZM_normalization4AT1_2 = 0;
end

p = 2;
a3 = 0.0;
if (strcmp(CZM_model_name, 'Linear') == 1)
    a2 = -0.5;
elseif (strcmp(CZM_model_name, 'Bilinear') == 1)
    a2 = 0.03687;
    a3 = 20.8343;
elseif (strcmp(CZM_model_name, 'Exponential') == 1)
    a2 = power(2.0, 5.0/3.0) - 3.0;
elseif (strcmp(CZM_model_name, 'Hyperbolic') == 1)
    a2 = power(2.0, 7.0/3.0) - 4.5;
elseif (strcmp(CZM_model_name, 'Concrete') == 1)
    a2 = 1.3868;
    a3 = 0.6567;
end

alphaFactor = 3.0 / 4.0; % sigma_f = 1, G = 0.5
alphaFactor = 3.0 / 8.0; % sigma_f = 2, G = 1
%alphaFactor is the reduction factor relative to 8/3/pi value provided in
%Phu
%alphaFactor > 0.75 -> invalid - sigma_f < 1 !
%alphaFactor = 0.75 -> sigma_f = 1 -> Gc = 0.5
%alphaFactor < 0.75 -> sigmaF = 0.75 / alphaFactor > 1
%alphaFactor = 3/8  -> sigmaF = 3/4/(3/8) = 2.0 -> Gc = 1.0
if (bPrime < 0)
    if (omegaCZM)
        bPrime = alphaFactor * 8/3/pi; % for CZM
    else
        if (CZM_normalization4AT1_2 == 1) % length normalization = l_coh
            if (xi == 0) % AT2
                bPrime = 27.0/256.0;
            elseif (xi == 1) % AT1
                bPrime = 3.0/8.0;
            end
        else
            bPrime = 1; % length normalization = b
        end
    end
end   



a23 = a2 * a3;
a1 = 4.0 / pi / bPrime;
Dmax = 1 - delD;

eps_p_i = 1;
if (xi == 0) % AT2
    c_alpha = 2.0;
    eps_p_i = 0;
elseif (xi == 1) % AT1
    c_alpha = 8.0/3.0;
    if (CZM_normalization4AT1_2 == 1)
        eps_p_i = sqrt(3.0/8.0);
    end
elseif (xi == 2) % CZM
    c_alpha = pi;
end

fomxi = 4.0 * (1.0- xi); 
D_vec = 0:Dmax/10000:Dmax;
sz = length(D_vec);
for i = 1:sz
    Dsn = D_vec(i);
    omD = 1.0 - Dsn;
    omega = omD * omD;
    omegaPrime = -2.0 * omD;
    if (omegaCZM)
        Ap = power(omD, p - 1);
        A = omD * Ap;
        Ap = -p * Ap;
        B = Dsn * (1.0 + a2 * Dsn + a23 * Dsn * Dsn);
        Bp = 1.0 + 2.0 * a2 * Dsn + 3.0 * a23 * Dsn * Dsn; 
        tmp = (A + a1 * B);
        omega = A / tmp;
        omegaPrime = -a1 * (A * Bp - Ap * B) / tmp / tmp;
    end
    eps = sqrt((2.0 * xi + fomxi * Dsn) / c_alpha / bPrime / -omegaPrime);
    sig = eps * omega;
    epsilon_p_vec(i) = eps;
    sigma_p_vec(i) = sig;
    omegaD_vec(i) = omega;
end

[sigmap_Max, i] = max(sigma_p_vec);
sigmap_Max_eps = epsilon_p_vec(i);
eps_f = epsilon_p_vec(sz);
phi = trapz(epsilon_p_vec, sigma_p_vec) + 0.5 * eps_p_i * eps_p_i;
phi_unloading = trapz(epsilon_p_vec(i:sz), sigma_p_vec(i:sz));
phi_loading = phi - phi_unloading;
brittleness_phi = phi_loading / phi;
brittleness_strain = sigmap_Max_eps / eps_f;

if (xi > 0)
    epsilon_p_vec = [[0] epsilon_p_vec];
    D_vec = [[0] D_vec];
    sigma_p_vec = [[0] sigma_p_vec];
    omegaD_vec = [[1] omegaD_vec];
end

delD = D_vec(2) - D_vec(1);
dDinv = 1.0 / delD;
dDiv2 = dDinv * dDinv;
szz = length(D_vec);
D_eps_dD = 0 * epsilon_p_vec;
D2_eps_dD2 = D_eps_dD;
D_eps_dD(1) = (epsilon_p_vec(2) - epsilon_p_vec(1)) * dDinv;
D_eps_dD(szz) = (epsilon_p_vec(szz) - epsilon_p_vec(szz - 1)) * dDinv;
for i = 2:szz-1
    D_eps_dD(i) = (epsilon_p_vec(i + 1) - epsilon_p_vec(i)) * dDinv;
    D2_eps_dD2(i) = (epsilon_p_vec(i + 1) + epsilon_p_vec(i - 1) - 2.0 * epsilon_p_vec(i)) * dDiv2;
end
D2_eps_dD2(1) = D2_eps_dD2(3);
D2_eps_dD2(2) = D2_eps_dD2(3);
D2_eps_dD2(szz) = D2_eps_dD2(szz - 2);
D2_eps_dD2(szz - 1) = D2_eps_dD2(szz - 2);

Dp_vec = 1.0 ./ D_eps_dD;
Dpp_vec = -D2_eps_dD2 .* Dp_vec.* Dp_vec.* Dp_vec;
Dpp_vec(1) = Dpp_vec(3);
Dpp_vec(2) = Dpp_vec(3);
Dpp_vec(szz) = Dpp_vec(szz - 2);
Dpp_vec(szz - 1) = Dpp_vec(szz - 2);


if (0)
    figure(1);
    plot(epsilon_p_vec, sigma_p_vec);
    figure(2);
    plot(epsilon_p_vec, D_vec);
    a = 12;
end