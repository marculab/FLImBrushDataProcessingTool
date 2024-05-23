idx1 = 6095;
% idx2 = 5809;
idx2 = 6606;
fit1 = get(Ch1DataObj,'fit',idx1);
fit2 = get(Ch1DataObj,'fit',idx2);
data = get(Ch1DataObj,'wf_aligned');

figure(1)
plot(data(:,idx1))
hold on
plot(fit1)
title(sprintf('Chi2: %.2f, Gain: %.2f',Ch1DataObj.stat_test.chi2.stat(idx1),Ch1DataObj.gain(idx1)))
hold off

figure(2)
plot(data(:,idx2))
hold on
plot(fit2)
title(sprintf('Chi2: %.2f, Gain: %.2f',Ch1DataObj.stat_test.chi2.stat(idx2),Ch1DataObj.gain(idx2)))
hold off

figure(3)
plot(Ch1DataObj.averagedData(:,idx1))
hold on
plot(Ch1DataObj.averagedData(:,idx2));
grid on
hold off
Ch1DataObj.gain(idx1)
Ch1DataObj.gain(idx2)

figure(4)
tiledlayout(2,1)
nexttile
plot(Ch1DataObj.stat_test.chi2.stat)
ylim([0 20])
title('Chi2 trace')
nexttile    
histogram(Ch1DataObj.stat_test.chi2.stat,[0:0.1:10])
title('Chi2 Histogram')