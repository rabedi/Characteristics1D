function [eps_vec, sig_vec, epsM, sigM, epsF, phiF] = ReadGiangData_NonDimIO(l10a, ldf, isHyper)
if (nargin < 1)
    l10a = 2.5;
end
if (nargin < 2)
    ldf = 0;
end
if (nargin < 3)
    isHyper = 1;
end
l10aGiang = l10a + 3;
[eps_vec, sig_vec, epsM, sigM, epsF, phiF] = ReadGiangData(l10aGiang, ldf, isHyper);
eps_vec = 0.1 * eps_vec;
sig_vec = 0.1 * sig_vec;
epsM = 0.1 * epsM;
sigM = 0.1 * sigM;
epsF = 0.1 * epsF;
phiF = 0.01 * phiF;

plot(eps_vec, sig_vec);
