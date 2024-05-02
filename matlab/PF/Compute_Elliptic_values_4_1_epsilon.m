function [vecVals, fullDamage] = Compute_Elliptic_values_4_1_epsilon(xi, eps, bBar)
% computes the vector of sigma, D, Dp Dpp, omega, and [epsilon = eps], for
% a given epsilon

% xi = 0 -> AT2, xi = 1 -> AT1, xi = 2, CZM-W
fullDamage = 0;
Dp = 0.0;
Dpp = 0.0;
if (xi == 2)
    epsFinal = 2.0 / pi / bBar;
    if (eps - epsFinal > -1e-6 * epsFinal)
        fullDamage = 1;
        sigma = 0.0;
        D = 1.0;
        omega = 0.0;
    elseif (eps < 1.0)
        sigma = eps;
        Dpp = 0.0;
        omega = 1.0;
    else
        sigma = (epsFinal - eps) / (epsFinal - 1.0);
        omega = sigma / eps;
        r = 1.0 / omega;
        a1 = 4.0 / pi / bBar;
        A = (r - 1 + 0.5 * a1);
        B = (2.0 - a1 - 2.0 * r);
        C = (r - 1.0);
        disc = B * B - 4.0 * A * C;
        if (disc < 0)
            disc
            fprintf(1, 'disc is negative\n');
            pause;
        end
        discr = sqrt(disc);
        fact = 1.0 / 2.0 / A;
        D1 = (-B - discr) * fact;
        D2 = (-B + discr) * fact;
        D1good = ((D1 >= -0.0001) && (D1 <= 1.00001));
        D2good = ((D2 >= -0.0001) && (D2 <= 1.00001));
        if (D1good)
            if (D2good)
                D = max(D1, D2);
            else
                D = D1;
            end
        else
            if (D2good)
                D = D2;
            else
                D1good
                D2good
                D1
                D2
                A
                B
                C
                eps
                epsFinal
                sigma
                omega
                fprintf(1, 'none of the Ds is acceptable');
                pause;
            end
        end
        D = min(max(D, 0.0), 1.0);
    end
else
    if (xi == 1) % AT1
        epsi2 = 3.0 / 8.0 / bBar;
        epsi = sqrt(epsi2);
        if (eps < epsi)
            D = 0.0;
            sigma = eps;
            omega = 1.0;
        else
            eps2 = eps * eps;
            tmp = epsi2 / eps2;
            D = 1.0 - tmp;
            sigma = tmp * epsi2 / eps;
            Dp = -2.0 * tmp / eps;
            Dpp = -3.0 * Dp / eps;
        end
    else % xi == 0 % AT2
        beps2 = bBar * eps * eps;
        opbeps2Inv = 1.0 / (1.0 + beps2);
        opbeps2Inv2 = opbeps2Inv * opbeps2Inv;
        D = beps2 * opbeps2Inv;
        sigma = eps * opbeps2Inv * opbeps2Inv2;
        Dp = 2.0 * bBar * sigma;
        Dpp = 2.0 * bBar * (3.0 * beps2 - 1.0) * opbeps2Inv * opbeps2Inv2;
    end
    omega = (1.0 - D);
    omega = omega * omega;
end
vecVals = [sigma, D, Dp, Dpp, omega, eps];
