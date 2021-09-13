%show time window maters
%%
clear all
close all
clc

%% Laguerre functions
order = 12;
alpha = 0.916;

%% 12.5GS/s
figure('Position',[50 50 1600 900])
tiledlayout(4,1)
N = 680;
L=Laguerre(680,order,alpha);
tau = 4;
dt = 0.08;
T = 54.4; % ns
t = 0:dt:T-dt;
t=t';
D = exp(-t/tau);
nexttile
plot(D,'r.--')
hold on
plot(L)
lt = h_lifet(D,dt,'average')
title(sprintf('Sampling rate %1f GS/s, # of data points: %d, averge lifetime: %.6f',1/dt, length(t),lt))
set(gca,'FontSize',12)
% 10GS/s
tau = 4;
dt = 0.1;
t = 0:dt:T-dt;
t=t';
D = exp(-t/tau);
nexttile
plot(D,'r.--')
hold on
plot(L)
lt = h_lifet(D,dt,'average')
title(sprintf('Sampling rate %1f GS/s, # of data points: %d, averge lifetime: %.6f',1/dt, length(t),lt))
set(gca,'FontSize',12)
% 5GS/s
tau = 4;
dt = 0.2;
t = 0:dt:T-dt;
t=t';
D = exp(-t/tau);
nexttile
plot(D,'r.--')
hold on
plot(L)
lt = h_lifet(D,dt,'average')
title(sprintf('Sampling rate %1f GS/s, # of data points: %d, averge lifetime: %.6f',1/dt, length(t),lt))
set(gca,'FontSize',12)
% 2.5GS/s
tau = 4;
dt = 0.4;
t = 0:dt:T-dt;
t=t';
D = exp(-t/tau);
nexttile
plot(D,'r.--')
hold on
plot(L)
lt = h_lifet(D,dt,'average')
title(sprintf('Sampling rate %1f GS/s, # of data points: %d, averge lifetime: %.6f',1/dt, length(t),lt))
set(gca,'FontSize',12)

%%

figure('Position',[50 50 1600 900])
tiledlayout(4,1)

tau = 4;
dt = 0.08;
T = dt*680; % ns
t = 0:dt:T-dt;
t=t';
D = exp(-t/tau);
nexttile
plot(t,D,'r.--')
hold on
L=Laguerre(length(t),order,alpha);
plot(t,L)
lt = h_lifet(D,dt,'average')
title(sprintf('Sampling rate %1f GS/s, # of data points: %d, averge lifetime: %.6f',1/dt, length(t),lt))
set(gca,'FontSize',12)
% 10GS/s
tau = 4;
dt = 0.1;
T = dt*680; % ns
t = 0:dt:T-dt;
t=t';
D = exp(-t/tau);
nexttile
plot(t,D,'r.--')
hold on
L=Laguerre(length(t),order,alpha);
plot(t,L)
lt = h_lifet(D,dt,'average')
title(sprintf('Sampling rate %1f GS/s, # of data points: %d, averge lifetime: %.6f',1/dt, length(t),lt))
set(gca,'FontSize',12)
% 5GS/s
tau = 4;
dt = 0.2;
T = dt*680; % ns
t = 0:dt:T-dt;
t=t';
D = exp(-t/tau);
nexttile
plot(t,D,'r.--')
hold on
L=Laguerre(length(t),order,alpha);
plot(t,L)
lt = h_lifet(D,dt,'average')
title(sprintf('Sampling rate %1f GS/s, # of data points: %d, averge lifetime: %.6f',1/dt, length(t),lt))
set(gca,'FontSize',12)
% 2.5GS/s
tau = 4;
dt = 0.4;
T = dt*680; % ns
t = 0:dt:T-dt;
t=t';
D = exp(-t/tau);
nexttile
plot(t,D,'r.--')
hold on
L=Laguerre(length(t),order,alpha);
plot(t,L)
lt = h_lifet(D,dt,'average')
title(sprintf('Sampling rate %1f GS/s, # of data points: %d, averge lifetime: %.6f',1/dt, length(t),lt))
set(gca,'FontSize',12)