function cdf = Stat_CDF_Weibull(x, x0, alpha, bTerm)
% CDF = 1 - exp(- bTerm (x - x0)^alpha)
if (x < x0)
	cdf = 0.0;
else
    delx = (x - x0);
    expTerm = bTerm * power(delx, alpha);
    cdf = 1.0 - exp(-expTerm);
end