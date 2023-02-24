d = load('xPos_2_InterfaceRawFinalSln.txt', '-ascii');
Z = 4.0;
figure(1);
plot(d(:,6) - d(:,5), d(:,2));
print('-dpng', 'xPos_0_InterfaceRawFinalSln_delu_s.png');
figure(2);
plot(d(:,1), d(:,2));
print('-dpng', 'xPos_0_InterfaceRawFinalSln_t_s.png');
figure(3);
plot(d(:,1), d(:,6) - d(:,5));
print('-dpng', 'xPos_0_InterfaceRawFinalSln_t_delu.png');

figure(4);
plot(d(:,1), d(:,2) - Z * d(:,3));
print('-dpng', 'xPos_0_InterfaceRawFinalSln_t_w.png');

figure(5);
t = d(:,1);
sz = length(t);
delt = t(2:sz) - t(1:sz - 1);
plot(1:sz - 1, log(delt))
print('-dpng', 'xPos_0_InterfaceRawFinalSln_log_delt.png');

