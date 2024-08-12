rho = 1;
E = 1;
G = 1;
b = 0.01;
epsilonDot = 1.0;
epsilonDotBar = epsilonDot;
bbar = b;
S = 10^13;
lnS = log(S);

eta = bbar * epsilonDotBar
%eta = epsilonDot * sqrt(b * b * b * rho / G)
%eta = 0.01

omdS = 1.0 / (1.0 + eta * lnS + 0.5 * eta * eta * lnS * log(S - 1));
omdS = omdS * omdS;

dStar = 1.0 - omdS
tStarBar = 1 / epsilonDotBar / sqrt(omdS)
saveBar = pi * bbar * (2.0 - log(S - 1) * eta)

