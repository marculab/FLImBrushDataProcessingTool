% code to test Laguerre function with different sampling rate
%% 
clear
close all
clc

%%
window = 680*0.08; %ns
dtV4 = 0.08;
dtFB = 0.1;

NV4 = window/dtV4;
NFB = window/dtFB;
tV4 = (0:NV4-1)*dtV4;
tFB = (0:NFB-1)*dtFB;

order = 12;

alphaFB = alpha_up(NFB,order);

LaguerreFB = Laguerre(NFB,order,alphaFB);
LaguerreV4 = interp1(tFB,LaguerreFB,tV4,'pchip');
%%
K=11;
figure
plot(tFB,LaguerreFB(:,K),'r.')
hold on
plot(tV4,LaguerreV4(:,K),'go')
%%

tV4 = (1:NV4)*dtV4;
tFB = (1:NFB)*dtFB;
%% plot Laguerre functions
K = 5;
figure
plot(tV4,LaguerreV4(:,K))
hold on
plot(tFB,LaguerreFB(:,K))
hold off
