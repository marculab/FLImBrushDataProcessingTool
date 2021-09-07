% code to plot laguerre functions with same time window but differenrt samplong rate
clear all
close all
clc
%%
dtV = [0.05;0.1;0.2;0.4];
tWindow =154.4; % ns
figure

for i = 1:length(dtV)
   dt = dtV(i);
   t = 0:dt:tWindow;
   K = 12;
   alpha=alpha_up(length(t),K,0.9);
   L=Laguerre(length(t),K,alpha);
   L = L./L(1,:);
   plot(t,L(:,2))
   hold on
end
hold off