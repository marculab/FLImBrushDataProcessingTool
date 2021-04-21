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
%% plot first 8 LBFs
figure('Position',[200 200 800 600])
plot(laguerreFuncsLong,'LineWidth',1.3)
title(sprintf('First 8 Laguerre basis functions, alpha = %.3f, k = %d',alphaUpperLim,L))
box off
grid on
legend('Order = 1','Order = 2','Order = 3','Order = 4','Order = 5','Order = 6','Order = 7','Order = 8')
xlabel('k')
ylabel('Amplitude (a.u.)')
set(gca,'LineWidth',1.3)
set(gca,'FontSize',15)
saveas(gcf,'Laguerre Functions','png')

%%
alpha = linspace(0.6,alphaUpperLim,4);
figure('Position',[200 200 800 1000])
tiledlayout(4,1)
for i = 1:length(alpha)
    nexttile
    laguerreFuncs = Laguerre(L,8,alpha(i));
    plot(laguerreFuncs(:,7),'LineWidth',1.3)
    ylim([-0.4 0.4])
    set(gca,'YTick',[-0.4 -0.2 0 0.2 0.4])
    if i~=4
    xticklabels({})
    end
    title(['l = 7th Laguerre function with $\alpha = $ ' sprintf('%.3f', alpha(i))],'Interpreter','latex')
    box off
    grid on
    set(gca,'LineWidth',1.3)
    set(gca,'FontSize',15)
end
box off
grid on
% legend('Order = 1','Order = 2','Order = 3','Order = 4','Order = 5','Order = 6','Order = 7','Order = 8')
xlabel('k')
saveas(gcf,'Laguerre Functions varying alpha','png')


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
