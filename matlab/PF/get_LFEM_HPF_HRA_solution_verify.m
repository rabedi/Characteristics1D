function [sigmap_Max_eps, sigmap_Max, phi_loading, eps_f, phi] = get_LFEM_HPF_HRA_solution_verify(model_s, CZM_normalization4AT1_2)

if (nargin < 1)
    model_s = 'AT1';
%    model_s = 'AT2';
end
if (nargin < 2)
    % 0: no normalization: lScale = b
    % 1: quasi-static strength = 1 
    CZM_normalization4AT1_2 = 1;
%    CZM_normalization4AT1_2 = 0;
end

[F, A, beta] = get_LFEM_PH_F_A_beta(model_s, CZM_normalization4AT1_2);
eps_m = 1.087417092379920/F^0.25
sig_m = 0.854307960622323/F^0.25
phi_m = 0.547542374099367/F^0.5
eps_f = 2.003147359426884/F^0.25
phi_f = 0.975198395496296/F^0.5
