function [xs_tri, pdf_gen_tri, cdf_gen_tri, xs_Weibull, pdf_Weibull, cdf_Weibull] = GetPW_Average_Distributions(alpha, delta, N)
if (nargin < 1)
    alpha = 2.0;
end
if (nargin < 2)
    delta = 0.5;
end
if (nargin < 3)
    N = 100;
end
x0 = 1 - delta;
xMax = 1.0 + delta;
xs_tri = x0:delta/100:xMax;
len = length(xs_tri);
for i = 1:len
    x = xs_tri(i);
    cdf_gen_tri(i) = Stat_CDF_GenTriangular(x, alpha, delta);
    pdf_gen_tri(i) = Stat_PDF_GenTriangular(x, alpha, delta);
end

xMax = 1.0 * xMax;
xs_Weibull = x0:delta/100:xMax;
len = length(xs_Weibull);
bTerm = N * 0.5 / power(delta, alpha);
for i = 1:len
    x = xs_Weibull(i);
    cdf_Weibull(i) = Stat_CDF_Weibull(x, x0, alpha, bTerm);
    pdf_Weibull(i) = Stat_PDF_Weibull(x, x0, alpha, bTerm);
end
