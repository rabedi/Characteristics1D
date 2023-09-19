function pdf = Stat_PDF_Weibull(x, x0, alpha, bTerm)
% PDF = CDF', CDF = 1 - exp(- bTerm (x - x0)^alpha) ->
% PDF = bTerm (x - x0)^(alpha - 1) * exp(- bTerm (x - x0)^alpha)
if (x < x0)
	pdf = 0.0;
else
    delx = (x - x0);
    expTerm = bTerm * power(delx, alpha - 1);
    pdf = alpha * expTerm * exp(-expTerm * delx);
end