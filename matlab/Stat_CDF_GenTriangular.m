function cdf = Stat_CDF_GenTriangular(x, alpha, delta)

if (x < (1.0 - delta))
	cdf = 0.0;
elseif (x > (1.0 + delta))
	cdf = 1.0;
else
    inv_2_delta_pow_alpha = 0.5 / power(delta, alpha);
    if (x < 1.0)
    	cdf = inv_2_delta_pow_alpha * power(x - (1.0 - delta), alpha);
    else
        cdf = 1.0 - inv_2_delta_pow_alpha * power((1.0 + delta) - x, alpha);
    end
end