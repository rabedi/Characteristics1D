function [yEpsp_vec, sigmap_vec, D_vec, Dp_vec,  Dpp_vec, omegaD_vec, sigmap_Max, sigmap_Max_y, yEpsp_f, ...
    phi, phi_unloading, phi_loading, brittleness_phi, brittleness_strain] = ...
    ComputePF(useAnalyticalSln, isHyper, xi, omegaCZM, bPrime, CZM_model_name, CZM_normalization4AT1_2)
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
    useAnalyticalSln = 0;
end
if nargin < 2
    isHyper = 0;
end
if nargin < 3
    xi = 2;
end
if nargin < 4
    if (xi == 2)
        omegaCZM = 1;
    else
        omegaCZM = 0;
    end
end

if nargin < 5
    bPrime = -1;
end

if nargin < 6
    CZM_model_name = 'Linear';
end

if nargin < 7
    CZM_normalization4AT1_2 = 1;
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
        if (xi == 0) % AT2
            bPrime = 27.0/256.0;
        elseif (xi == 1) % AT1
            bPrime = 3.0/8.0;
        end
    end
end   


if (omegaCZM)
    if (isHyper)
        beta = power(1.0/3.0/bPrime^2, 1.0/3.0);
    else
        beta = power(1.0/3.0/bPrime, 1.0/3.0);
    end
else
    if (xi == 1) % AT1
        beta = power(2.0/3.0, 2.0/3.0);
    elseif (xi == 0) % AT2
        beta = power(1.0/3.0, 1.0/3.0);
    end
    if (CZM_normalization4AT1_2)
        beta = beta * power(bPrime, 1.0 / 3.0);
    end
end

if (isHyper)
    delD = 0;
end

a23 = a2 * a3;
a1 = 4.0 / pi / bPrime;
Dmax = 1 - delD;
dely2 = dely * dely;
tbeta3 = 3.0 * beta * beta * beta;

if (omegaCZM)
    useAnalyticalSln = 0;
end

if (isHyper)
    Ds = 0.0;
    y = 0.0;
    Dds = 0.0;
    Ddds = 0.0;
    cntr = 0;
    while (Ds < Dmax)
        cntr = cntr + 1;
        yEpsp_vec(cntr) = y;
        D_vec(cntr) = Ds;
        Dp_vec(cntr) = Dds;
        Dpp_vec(cntr) = Ddds;
        yNew = y + 0.5 * dely;
        Dsn = Ds + dely * Dds + 0.5 * dely2 * Ddds; 
        omD = 1.0 - Dsn;
        omegaPrime = omD;
        if (omegaCZM)
            Ap = power(omD, p - 1);
            A = omD * Ap;
            Ap = -p * Ap;
            B = Dsn * (1.0 + a2 * Dsn + a23 * Dsn * Dsn);
            Bp = 1.0 + 2.0 * a2 * Dsn + 3.0 * a23 * Dsn * Dsn; 
            tmp = (A + a1 * B);
            omegaPrime = (A * Bp - Ap * B) / tmp / tmp;
        end
        % DddsNew = yNew; % * yNew; % * tbeta3 * omegaPrime;
        % DddsNew = omD; % * yNew; % * tbeta3 * omegaPrime;
        DddsNew = yNew * yNew * tbeta3 * omegaPrime;
        DdsNew = Dds + ((1 - delta) * Ddds + delta * DddsNew) * dely;
        Ds = Dsn;
        Dds = DdsNew;
        Ddds = DddsNew;
        y = y + dely;
    end
else
    Ds = 0.0;
    y = 0.0;
    Dds = 0.0;
    Ddds = 0.0;
    cntr = 0;
    while (Ds < Dmax)
        cntr = cntr + 1;
        yEpsp_vec(cntr) = y;
        D_vec(cntr) = Ds;
        Dp_vec(cntr) = Dds;
        Dpp_vec(cntr) = Ddds;
        yNew = y + 0.5 * dely;
        Dsn = Ds; % + dely * Dds + 0.5 * dely2 * Ddds; 
        omD = 1.0 - Dsn;
        omegaPrime = omD;
        if (omegaCZM)
            Ap = power(omD, p - 1);
            A = omD * Ap;
            Ap = -p * Ap;
            B = Dsn * (1.0 + a2 * Dsn + a23 * Dsn * Dsn);
            Bp = 1.0 + 2.0 * a2 * Dsn + 3.0 * a23 * Dsn * Dsn; 
            tmp = (A + a1 * B);
            omegaPrime = (A * Bp - Ap * B) / tmp / tmp;
        end
        DsNew = Ds + omegaPrime * tbeta3 * yNew * yNew * dely;
        Ds = DsNew;
        yNew = y + dely;
        DdsNew = omegaPrime * tbeta3 * yNew * yNew;
        Ddds = (DdsNew - Dds) / dely;
        Dds = DdsNew;
        y = yNew;
    end
end
%plot(yEpsp_vec, D_vec);
sz = length(yEpsp_vec);
sigmap_vec = 0 * D_vec;
omegaD_vec = sigmap_vec; 
for i = 1:sz
    y = yEpsp_vec(i);
    Dsn = D_vec(i);
    omD = 1.0 - Dsn;
    omega = omD * omD;
    if (omegaCZM)
        A = power(omD, p);
        B = Dsn * (1.0 + a2 * Dsn + a23 * Dsn * Dsn);
        tmp = (A + a1 * B);
        omega = A / tmp;
    end
    omegaD_vec(i) = omega;
    sigmap = omega * y;
    sigmap_vec(i) = sigmap;
end
[sigmap_Max, i] = max(sigmap_vec);
sigmap_Max_y = yEpsp_vec(i);
yEpsp_f = yEpsp_vec(sz);
phi = trapz(yEpsp_vec, sigmap_vec);
phi_unloading = trapz(yEpsp_vec(i:sz), sigmap_vec(i:sz));
phi_loading = phi - phi_unloading;
brittleness_phi = phi_loading / phi;
brittleness_strain = sigmap_Max_y / yEpsp_f;
if (0)
    figure(1);
    plot(yEpsp_vec, sigmap_vec);
    figure(2);
    plot(yEpsp_vec, D_vec);
    a = 12;
end