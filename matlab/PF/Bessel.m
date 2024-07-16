nu = 0.25;
format long;
zroot = fzero(@(z) besselj(-nu, z), 2)
besselj(-nu, zroot)

eps = 0.0000000000000001;
x = eps:eps:0.000000000001;
mod = power(x, nu);
% 1st
%y = besselj(nu,x);
% 2nd
y = mod .* bessely(nu,x);
% m 1st
%y = besseli(nu,x);
% m 2nd
%y = mod .* besselk(nu,x);
plot(x, y);