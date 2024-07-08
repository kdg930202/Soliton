%Turn on resonant pumping

[time, s08] = pth_change(1.2, 0.8);
[time, s1] = pth_change(1.2, 1);
[time, s2] = pth_change(1.2, 2);
[time, s4] = pth_change(1.2, 4);
save("homo_nonresonant.mat")



% time = 0:0.1:500;

figure()
plot(time, s08, 'LineWidth', 2)
hold on
plot(time, s1, 'LineWidth', 2)
plot(time, s2, 'LineWidth', 2)
plot(time, s4, 'LineWidth', 2)
axis tight
legend({'0.8Pth','1Pth','2Pth','4Pth'},'FontSize',14)
xlabel('Time(ps)')