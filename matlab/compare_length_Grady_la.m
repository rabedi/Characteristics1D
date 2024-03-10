x = -3:0.1:4;
y1 = 1.0 / 3.0 * log10(24.0) - 2.0 / 3.0 * x;
y2 = -x;
plot(x, y1, '-b');
hold on;
plot(x, y2, '-r');
legend({'Grady', 'la'});