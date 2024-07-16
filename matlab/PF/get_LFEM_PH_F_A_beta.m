function [F, A, beta] = get_LFEM_PH_F_A_beta(model_s, CZM_normalization4AT1_2)

% F = (ca * bBar / 2) = 4 * A * A = 3 * beta^3
% A = sqrt(F)/2
% beta = power(F/3, 1/3)
if (nargin < 1)
    model_s = 'AT1';
end
if (nargin < 2)
    % 0: no normalization: lScale = b
    % 1: quasi-static strength = 1 
    % 2: A is made equal to 1: useful for HRA-HPF
    CZM_normalization4AT1_2 = 1;
end

if (CZM_normalization4AT1_2 == 2)
    F = 4.0;
else
    if (strcmp(model_s, 'AT1') == 1)
        if (CZM_normalization4AT1_2 == 0)
            F = 4.0 / 3.0;
        else
            F = 0.5;
        end
    else
        if (CZM_normalization4AT1_2 == 0)
            F = 1.0;
        else
            F = 27.0 / 256.0;
        end
    end
end

A = 0.5 * sqrt(F);
beta = power(F/3.0, 1/3.0);
