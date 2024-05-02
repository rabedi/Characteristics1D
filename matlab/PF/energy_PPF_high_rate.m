isAT1 = 1;
scaleBbarAsCZM = 1;
% D = 1 - tol is considered full damage
tol = 0.001;

if (isAT1)
    if (scaleBbarAsCZM)
        beta = power(1.0/6.0, 1.0/3.0);
    else
        beta = power(2.0/3.0, 2.0/3.0);
    end
else
    if (scaleBbarAsCZM)
        beta = power(3.0/16.0, 2.0/3.0);
    else
        beta = power(1.0/3.0, 1.0/3.0);
    end
end

tMax = 1.0 / 3.0;
yMax = power(tMax / 2.0, 1.0 / 3.0) / beta
sigmaMax = exp(-1.0 / 3.0) * yMax
z = 2.0 / 3.0;
factor = power(2.0, -2.0 / 3.0) / beta / beta / 3.0;

gi = gammainc(tMax, z) * gamma(z)
ene_0_to_max = factor * gi
gc = gamma(z)
ene_0_to_f = factor * gc
a = 10;

tol
tf = 2 * log(1/tol);
yf = power(tf / 2.0, 1.0 / 3.0) / beta
