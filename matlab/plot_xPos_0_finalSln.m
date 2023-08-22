%d = load('xPos_2_InterfaceRawFinalSln.txt', '-ascii');
addExact = 1;
a = 10.0;
d = load('../run_0xPos_0_InterfaceRawFinalSln.txt', '-ascii');
Z = 4.0;
Z = 1.0;
figure(1);
%delu = d(:,6) - d(:,5);
delu = d(:,7) - d(:,6);
plot(delu, d(:,2));
print('-dpng', 'xPos_0_InterfaceRawFinalSln_delu_s.png');
figure(2);
time = d(:,1);
plot(time, d(:,2));
print('-dpng', 'xPos_0_InterfaceRawFinalSln_t_s.png');
figure(3);
plot(time, delu);
if (addExact)
    deluExact = a * (0.5 * (exp(2.0 * time) - 1) - time);
    hold on;
    for j = 1:length(deluExact)
        if (deluExact(j) > 1)
            deluExact(j) = nan;
        end
    end
    plot(time, deluExact, '-r');
end
print('-dpng', 'xPos_0_InterfaceRawFinalSln_t_delu.png');


figure(4);
%plot(d(:,1), d(:,2) - Z * d(:,3));
plot(d(:,1), d(:,2) - Z * d(:,4));
print('-dpng', 'xPos_0_InterfaceRawFinalSln_t_w.png');

figure(5);
t = d(:,1);
sz = length(t);
delt = t(2:sz) - t(1:sz - 1);
plot(1:sz - 1, log(delt))
print('-dpng', 'xPos_0_InterfaceRawFinalSln_log_delt.png');

