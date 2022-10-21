figure;
t = [0:1/10:10]
plot(t,grandAvg{3}(:,4),'k','lineWidth',3)
hold on
plot(t,grandAvg{4}(:,4),'r','lineWidth',3)
lgd = legend('+4dB','-2dB');
lgd.FontSize = 20;
grid on;
xlabel('Time (sec)','FontSize',25)
xticks([0:1:15])
ylabel('Oxygenation (M.mm)','FontSize',25)
title('Left Frontal (Low Context)','FontSize',30)
xlim([0 15])