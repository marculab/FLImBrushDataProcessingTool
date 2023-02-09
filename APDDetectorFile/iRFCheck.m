i = 101;
step = 20;
%
figure
plot(irf(:,i)/max(irf(:,i)))
hold on
plot(irf(:,i+step)/max(irf(:,i+step)))
hold off
xlim([725 1300]/4)
irfV(i)
G = interp1(gainV,gain,irfV(i))
%%
figure
plot(irfUpSampled(:,i)/max(irfUpSampled(:,i)))
hold on
plot(circshift(irfUpSampled(:,i+step)/max(irfUpSampled(:,i+step)),2))
hold off
xlim([725 1300])

