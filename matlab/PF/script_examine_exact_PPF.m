xi = 1;
la_s = '2';
la = str2double(la_s);
aTimesf = power(10, la);
bPrime = -1;
CZM_normalization4AT1_2 = 1;
fastComputation = 0;
delD = -1;
if (xi == 1)
    model_s = 'AT1';
else
    model_s = 'AT2';
end

isHyper = 0;
isApproximate = 0;
l_cD2c_s = '0';
df_s = '1';
bTimesPiCZM_s = '1';
CZM_modelName = 'Linear';
tauModel = 'b2cd';

PFtmp = PF;
vals_default = {isHyper, isApproximate, model_s, la_s, l_cD2c_s, df_s, CZM_normalization4AT1_2, bTimesPiCZM_s, CZM_modelName, tauModel};
[PFtmp, vals_out] = PFtmp.Initialize_Stage1(vals_default);
PFtmp = PFtmp.Compute();


epsilon_p_vec{2} = PFtmp.epsilon_p;
sigma_p_vec{2} = PFtmp.sigma_p;
D_vec{2} = PFtmp.D_vec;
Dp_vec{2} = PFtmp.Dp_vec;
Dpp_vec{2} = PFtmp.Dpp_vec;
omegaD_vec{2} = PFtmp.omegaD_vec;
sigmap_Max(2) = PFtmp.sigmap_Max;
sigmap_Max_eps(2) = PFtmp.sigmap_Max_eps;
eps_f(2) = PFtmp.eps_f;
phi(2) = PFtmp.phi;
phi_unloading(2) = PFtmp.phi_unloading;
phi_loading(2) = PFtmp.phi_loading;
brittleness_phi(2) = PFtmp.brittleness_phi;
brittleness_strain(2) = PFtmp.brittleness_strain;

[epsilon_p_vec{1}, sigma_p_vec{1}, D_vec{1}, Dp_vec{1}, Dpp_vec{1}, omegaD_vec{1}, sigmap_Max(1), sigmap_Max_eps(1), eps_f(1), ...
    phi(1), phi_unloading(1), phi_loading(1), brittleness_phi(1), brittleness_strain(1)] = ...
    ComputePF_ParabolicAT12(xi, aTimesf, bPrime, CZM_normalization4AT1_2, fastComputation, delD);


clrs = {'r', 'b'};
figure(1);
for i = 1:2
    plot(epsilon_p_vec{i}, sigma_p_vec{i}, 'Color', clrs{i}, 'LineWidth', 2);
    hold on;
end

figure(2);
for i = 1:2
    plot(epsilon_p_vec{i}, D_vec{i}, 'Color', clrs{i}, 'LineWidth', 2);
    hold on;
end

sigmap_Max
phi
