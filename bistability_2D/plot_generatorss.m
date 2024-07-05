%Turn off resonant pumping
% s08 = pth_change(0.8);
% s1 = pth_change(1);
% s2 = pth_change(2);
% s4 = pth_change(4);
% save("tau_p = 15ps.mat")

time = 0:0.1:500;
figure()
plot(time, s08, 'LineWidth', 2)
hold on
plot(time, s1, 'LineWidth', 2)
plot(time, s2, 'LineWidth', 2)
plot(time, s4, 'LineWidth', 2)
axis tight
legend({'0.8Pth','1Pth','2Pth','4Pth'},'FontSize',14)
xlabel('Time(ps)')