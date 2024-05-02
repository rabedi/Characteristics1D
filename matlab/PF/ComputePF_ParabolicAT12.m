function [epsilon_p_vec, sigma_p_vec, D_vec, Dp_vec, Dpp_vec, omegaD_vec, sigmap_Max, sigmap_Max_eps, eps_f, ...
    phi, phi_unloading, phi_loading, brittleness_phi, brittleness_strain] = ...
    ComputePF_ParabolicAT12(xi, aTimesf, bPrime, CZM_normalization4AT1_2, fastComputation, delD)
%CZM_normalization4AT1_2: means that when we use AT1/AT2 models, we change
%the normalization from sigmaNormalization = sqrt(G/Eb) and
%lengthNormalization = b to sigmaNormalization = sigma_coh, and
%lengthNormalization = l_coh. This option = 1 is useful if we compare AT1
%and AT2 with CZM models



if nargin < 1
    xi = 1;
end
if nargin < 2
    % loading rate times f
    aTimesf = 1.0;
    aTimesf = 0.01;
end

if nargin < 3
    bPrime = -1;
end

if nargin < 4
    CZM_normalization4AT1_2 = 0;
%    CZM_normalization4AT1_2 = 1;
end

if nargin < 5
    fastComputation = -1;
end

if nargin < 6
    delD = -1;
end

if (delD < 0)
    delD = 1e-3;
    delD = 3e-2;
    delD = 6e-2;
end

% if 1 -> incrementally calculates the integral of exponential -> less
% accurate and blows up for low a
% else for each x it calculates the exponential integral
if (fastComputation == -1)
    if (aTimesf < 1e-2)
        fastComputation = 0;
    else
        fastComputation = 0; %1; 
    end
end

del_eps = 1e-4;
if (aTimesf > 10000)
    del_eps = max(2e-3, del_eps);
elseif (aTimesf > 1000)
    del_eps = max(3e-4, del_eps);
end

if (bPrime < 0)
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


Dmax = 1 - delD;

eps_p_i = 1;
if (xi == 0) % AT2
    c_alpha = 2.0;
    eps_p_i = 0;
elseif (xi == 1) % AT1
    c_alpha = 8.0/3.0;
    eps_p_i = sqrt(3.0/8.0/bPrime);
end
eps_p_i3 = eps_p_i * eps_p_i * eps_p_i;

a = aTimesf;
b = bPrime;
if (xi == 1)
    f = 4.0 * b / 9.0 / a;
end

fomxi = 4.0 * (1.0 - xi); 
epsv(1) = eps_p_i;
Dv(1) = 0.0;
D = 0.0;
cntr = 1;
s35 = sqrt(3.0 / 5.0);
pts(1) = 0.5 * (-s35 + 1);
pts(2) = 0.5;
pts(3) = 0.5 * (s35 + 1);
wts(1) = 5.0 / 18.0;
wts(3) = wts(1);
wts(2) = 4.0 / 9.0;

n = 3;
taInv = 1.0 / 2.0 / a;
IntVPrev = 0.0;
while (D < Dmax)
    epsp = epsv(cntr);
    Dp = Dv(cntr);
    epsn = epsp + del_eps;
    epsn3 = epsn * epsn * epsn;
    if (xi == 1)
        IntV = 0.0;
        if (fastComputation == 0)
            for j = 1:cntr
                epssbase = eps_p_i + (j - 1) * del_eps;
                for k = 1:n
                    w = wts(k) * del_eps;
                    epss = epssbase + pts(k) * del_eps;
                    epss3 = epss * epss * epss;
                    v = exp(-f * (epsn3 - epss3));
                    IntV = IntV + v * w;
                end
            end
        else
            j = cntr;
            IntVTmp = 0.0;
            epssbase = eps_p_i + (j - 1) * del_eps;
            for k = 1:n
                w = wts(k) * del_eps;
                epss = epssbase + pts(k) * del_eps;
                epss3 = epss * epss * epss;
                v = exp(f * epss3);
                IntVTmp = IntVTmp + v * w;
            end
            IntV = IntVPrev + IntVTmp * exp(-f * epsn3);
            IntVPrev = IntVPrev + IntVTmp;
        end
        D = 1.0 - exp(-f * (epsn3 - eps_p_i3)) - taInv * IntV;
    end
    cntr = cntr + 1;
    epsv(cntr) = epsn;
    Dv(cntr) = D;
end

epsilon_p_vec = epsv;
D_vec = Dv;
sz = length(D_vec);
for i = 1:sz
    D = D_vec(i);
    omd = 1 - D;
    omegaD_vec(i) = omd * omd;
    sigma_p_vec(i) = omd * omd * epsilon_p_vec(i);
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

if (0)
    figure(1);
    plot(epsilon_p_vec, sigma_p_vec);
    figure(2);
    plot(epsilon_p_vec, D_vec);
end
end