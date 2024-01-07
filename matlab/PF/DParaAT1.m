function D = DParaAT1(y)
beta3 = 4.0 / 9.0;
D = 1 - exp(-beta3 .* y .* y .* y); 
