function [eps_vec, sig_vec, epsM, sigM, epsF, phiF] = ReadGiangData(l10aGiang, ldf, isHyper)
if (nargin < 1)
    l10aGiang = 3.5;
end
if (nargin < 2)
    ldf = 1;
end
if (nargin < 3)
    isHyper = 1;
end
laGI = floor(l10aGiang + 1e-8);
laGR = l10aGiang - laGI;
if (laGR > 0)
    if (abs(laGR - 0.5) < 0.01)
        if (abs(l10aGiang + 0.5) < 1e-2)
            laGs = ['-0p5'];
        else
            laGs = [num2str(laGI), 'p5'];
        end
    else
        fprintf(1, 'Implement it here\n');
        pause
    end
else
    laGs = num2str(laGI);
end
ldfs = num2str(ldf);
if (isHyper == 1)
    DE = 'H';
elseif (isHyper == 0)
    DE = 'P';
else
    DE = 'E';
end

fileName = ['../../../Giang/bitraction_bar1D_', DE, 'PF_AT1_b0p01_a', laGs, '_df', ldfs, '_cl0_rf10/load_data.txt'];
dat = load(fileName);

%1 time
%2 loading rate
%3 sigma
%4 phi
%5 phi_r
%6 phi_d
%7 local bulk energy
%8 local crack energy
%9 phi
%10 10
ag = power(10.0, l10aGiang);
eps_vec = ag * dat(:, 1);
sig_vec = dat(:, 3);
[m, n] = size(dat);
%phiF = data(m, 6);
phiF = dat(m, 9);
epsF = eps_vec(m);
[sigM, i] = max(sig_vec);
epsM = eps_vec(i);
sigTol = 0.001 * sigM;
is = find(sig_vec > sigTol);
iss = length(is);
if (iss > 0)
    i = is(iss);
    epsF = eps_vec(i);
end
