function [slope, intercept, rsq, regs] = LinearRegression(x, y)

p = polyfit(x, y, 1);
slope = p(1);
intercept = p(2);
% Calculate R-squared value
yfit = polyval(p, x);
yresid = y - yfit;
SSresid = sum(yresid.^2);
SStotal = (length(y)-1) * var(y);
rsq = 1 - SSresid/SStotal;
regs = [slope, intercept, rsq];
