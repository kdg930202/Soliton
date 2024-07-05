close all
Pr_l = time_l*Inten_Pr;
Pr_r = time_r*Inten_Pr;


figure()
plot(Pr_l, l)
hold on
plot(Pr_r, r, 'r')
xlabel('Pr','FontSize',20)
xline(1.2,'--')
set(gca,'FontSize',20)