L = 2e-3;
%E = 680e9;
E = 610e9;
rho = 3900;
G = 83.13;
sigmac = 1e9;
h = L/2^14
b = 200e-9


lc = E * G / sigmac / sigmac
b2lc = b / lc
b2h = b / h
h2lc = h / lc
lc2h = lc / h
epsbar_f =  2 / pi / b2lc
cu = sqrt(E/rho)
epsDotc = sigmac^3 * cu / G / E / E
tauc = lc / cu
phic = G / lc
epsdTrans = 1e8 / epsDotc



%% Marigo's work
epsDot = 1;
b = 0.01;
E = 100;
rho = 100;
eta = epsDot * b^1.5 * sqrt(rho)
%eta = epsDot * b
S = 10^13;
lS = log(S);
lS1 = log(S - 1);
tmp = (1 + eta * lS + eta * eta * lS1);
ds = 1 - 1 / tmp / tmp
ts = 1 / epsDot / (1 - ds)
save = pi * b * (2 - lS1 * eta)