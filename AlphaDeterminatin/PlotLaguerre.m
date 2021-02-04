% code to plot Laguerre functinos of different truncation length and alpha
% value

%%
clear all
close all
clc

%% length = 1000;
O = 12;
L = 1000;
alphaUpperLim=alpha_up(L,O,[],[]);
laguerreFuncsLong = Laguerre(L,O,alphaUpperLim);

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
