PFtmp = PF;
PFtmp.omegaCZM = -1;
PFtmp.xi = 2;
bTimesPiCZM = 1.0; % for CZM-Wu value of 2 is max (sharp drop), 1 is good, 0.5 is two times slower decay, etc.
PFtmp.bTimesPiCZM = bTimesPiCZM;

PFtmp.ap = 100;
PFtmp.isHyper = 1;
% PFtmp.del_eps = -0.001;
PFtmp.del_eps = 0.001;
% asymptoticMode = 
%                       2  ->    high loading rate approximation AND only
%                       the asymptotic part
%                       1  ->    high loading rate approximation
%                       -1 ->   elliptic limit
%                       0  ->   compute the exact solution
PFtmp.asymptoticMode = 0;
plotResults = 0;
PFtmp.Compute(plotResults);