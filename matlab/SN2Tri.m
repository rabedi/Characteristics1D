function triVal = SN2Tri(xstandardNormalValue, dd2)
if nargin < 2
    dd2 = 0.5;
end

minV = 1 - dd2;
maxV = 1 + dd2;
inv_sqrt2 = 1.0 / sqrt(2.0);
p = 0.5 * (1.0 + erf(xstandardNormalValue * inv_sqrt2));
if (p < 0.5)
    triVal = 1.0 + dd2 * (-1.0 + sqrt(2.0 * p));
else
    triVal = 1 + dd2 * (1.0 - sqrt(2.0 * (1.0 - p)));
end
