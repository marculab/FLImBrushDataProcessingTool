% code to plot Laguerre functinos of different truncation length and alpha
% value

%%
clear all
close all
clc

%% length = 1000;
O = 8;
L = 500;
alphaUpperLim=alpha_up(L,O,[],[]);
laguerreFuncsLong = Laguerre(L,O,alphaUpperLim);
%%
figure('Position',[200 200 800 600])
plot(laguerreFuncsLong,'LineWidth',1.3)
title(sprintf('First 8 Laguerre basis functions, alpha = %.3f, k = %d',alphaUpperLim,L))
box off
grid on
xlabel('k')
ylabel('Amplitude (a.u.)')
set(gca,'LineWidth',1.3)
set(gca,'FontSize',15)
saveas(gcf,'Laguerre Functions','png')

%%
laguerreFuncsAlphaHalf = Laguerre(L,O,alphaUpperLim*0.8);
idx = 12;
figure
plot(laguerreFuncsLong(:,idx)/max(laguerreFuncsLong(:,idx)));
hold on
plot(laguerreFuncsAlphaHalf(:,idx)/max(laguerreFuncsAlphaHalf(:,idx)));
hold off
legend('Original','Half Alpha')


%% length by 2
alphaUpperLim2=alpha_up(L/2,O,[],[]);
laguerreFuncsShort = Laguerre(L/2,O,alphaUpperLim);
