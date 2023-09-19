function pdf = Stat_PDF_GenTriangular(x, alpha, delta)

if (x < (1.0 - delta))
	pdf = 0.0;
elseif (x > (1.0 + delta))
	pdf = 1.0;
else
    inv_2_delta_pow_alpha = 0.5 / power(delta, alpha);
    if (abs(alpha - 1.0) < 1e-3)
        pdf = inv_2_delta_pow_alpha;
    else
        if (x < 1.0)
	        pdf = alpha * inv_2_delta_pow_alpha * power(x - (1.0 - delta), alpha - 1.0);
        else
            pdf = alpha * inv_2_delta_pow_alpha * power((1.0 + delta) - x, alpha - 1.0);
        end
    end
end