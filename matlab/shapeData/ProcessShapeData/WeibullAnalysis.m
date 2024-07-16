function [m, lambda, N, R2, p, X, Y, xreg, yreg] = WeibullAnalysis(xs, xmin)
if (nargin  < 1)
    xs = wblrnd(2, 1, 100, 1); % Generate 100 random numbers from a Weibull distribution with shape 2 and scale 1
end
if (nargin < 2)
    xmin = 0.0;
end

xs = sort(xs) - xmin;
n = length(xs);
ys = 0.5/n:1/n:1;
ys = ys';
inds = find(xs <= 0);
linds = length(inds);
if (linds > 0)
    lastInd = inds(linds);
    xs = xs(lastInd + 1:n);
    ys = ys(lastInd + 1:n);
end
X = log(xs);
Y = log(-log((1 - ys)));


% Perform linear regression
p = polyfit(X, Y, 1); % p(1) is the slope, p(2) is the intercept

% Calculate the fitted values
yfit = polyval(p, X);

% Calculate the residuals
residuals = Y - yfit;

% Calculate the total sum of squares
SStot = sum((Y - mean(Y)).^2);

% Calculate the residual sum of squares
SSres = sum(residuals.^2);

% Calculate R-squared
R2 = 1 - (SSres / SStot);
xreg = [min(X), max(X)];
yreg = p(1) * xreg + p(2);

% p(1) is the slope, p(2) is the intercept
m = p(1);
C = p(2);
N = exp(C); % number N in Weibull distribution
lambda = power(N, -1.0/m);
