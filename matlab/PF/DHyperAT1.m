function D = DHyperAT1(y)
fprint(1, 'This formula is wrong!');
pause
pause

beta = power(2.0/3.0, 2.0/3.0);
beta3 = power(beta, 3.0);
D = y - gamma(1.0 / 3.0) * gammainc(beta3 * y.^3,1.0/3.0) / power(12, 1.0/3.0);