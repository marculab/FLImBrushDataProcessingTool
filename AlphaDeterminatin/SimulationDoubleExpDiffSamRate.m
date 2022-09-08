% code to run simulated deconvolution to determine the best alpha value

%%
clear 
close all
clc
addpath(genpath('C:\Users\Xiangnan\Documents\MyGitRepo\FLImBrushDataProcessingTool\Algorithms'))
%% get random lifetime
N = 5000;
% trueLT = rand(N,1);
tau1Low = 0.5;
tau1High = 4;
tau2 = 8;
tau1 = tau1Low+rand(N,1)*(tau1High-tau1Low);

k=9; % photon ratio

%% get FB decays at 10 GS/s
dtFB = 0.1;
tWindow = 54.4; % ns same for V4 nad V5
trueLTFormula = zeros(size(tau1));
t = 0:dtFB:tWindow-dtFB;
decayFB = zeros(round(tWindow/dtFB),N);
for i = 1:N
    M = [1 1;tau1(i) -k*tau2];
    A = M\[1;0];
    trueLTFormula(i) = (A(1)*tau1(i).^2+A(2)*tau2^2)./(A(1)*tau1(i)+A(2)*tau2);
    decayFB(:,i) = A(1)*exp(-t/tau1(i))+A(2)*exp(-t/tau2);
    
end
trueLTFB = h_lifet(decayFB,dtFB,'average');
binEdge = linspace(min(trueLTFB),max(trueLTFB),40);

figure
scatter(tau1,trueLTFormula-trueLTFB);
% title('histogram of randomly distributed lifetimes')
xlabel('tau1 (ns)')
ylabel('Lifetime difference')

figure
scatter(tau1,trueLTFormula);
title('Average lifetime vs tau1, fixed tau2 = 16 ns')
xlabel('tau1 (ns)')
ylabel('Lifetime')

figure
plot(decayFB(:,1:500:end))



%% get V4 decays at 12.5 GS/s
dtV4 = 0.08;
% tWindow = 54.4; % ns same for V4 nad V5
% trueLTC = zeros(size(tau1));
t = 0:dtV4:tWindow-dtV4;
decayV4 = zeros(round(tWindow/dtV4),N);
for i = 1:N
    M = [1 1;tau1(i) -k*tau2];
    A = M\[1;0];
%     trueLTC(i) = (A(1)*tau1(i).^2+A(2)*tau2^2)./(A(1)*tau1(i)+A(2)*tau2);
    decayV4(:,i) = A(1)*exp(-t/tau1(i))+A(2)*exp(-t/tau2);
    
end
trueLTV4 = h_lifet(decayV4,dtV4,'average');



%% load and truncate irf FB
irfStruc = load('..\APDDetectorFile\UCD_MC\M00549707_DCS.mat'); % APD irf
UpFactor = irfStruc.irfRawdt/dtFB;
irfFB = interp(irfStruc.irf(:,200),UpFactor);
[~,irfMaxIdx] = max(irfFB);
tLength = tWindow/dtFB;
startIdx = round(irfMaxIdx-0.2*tLength);
if startIdx <=0
    startIdx = 1;
end
try
irfTFB = irfFB(startIdx:startIdx+tLength-1);
catch
    irfTFB = irfFB(startIdx:end);
end
[~,MaxIdx] = max(irfTFB);
irfTFB = irfTFB./sum(irfTFB);
% irfTFB = circshift(irfTFB,266-MaxIdx);
figure
plot(irfTFB)
title('irf')

%% resample irf for V4
tFB = (0:tWindow/dtFB-1)*dtFB;
tFB = tFB';
tV4 = (0:tWindow/dtV4-1)*dtV4';
tV4 = tV4';

irfTV4 = interp1(tFB,irfTFB,tV4,'pchip');
irfTV4 = irfTV4./sum(irfTV4);
figure;plot(irfTV4)
%% convolusion FB
specFB = filter(irfTFB,1,decayFB);
specFB = specFB./max(specFB);
SNR = 55; %in dB, SNR = 20log10(Max/noise)
if SNR==0
    noise = zeros(size(specFB));
else
    noise = randn(size(specFB))*1/db2mag(SNR);
end
DC = 1:size(specFB,1);
DC = DC';
[~,MaxIdx] = max(specFB);
MaxIdx = mode(MaxIdx);
DC(DC<MaxIdx)=0;
DC(DC>=MaxIdx)=1;
DC = DC*0.0000;
specFB = specFB+noise+DC;
figure
plot(irfTFB)
hold on
plot(specFB(:,1:50:end))
hold off
title('Simulation data and irf')

%% deconvolution
LagOrder = 12;
alphaUpperLim=alpha_up(size(specFB,1),LagOrder,[],[])
% alphaUpperLim=0.916;

numOfAlpha = 1;
% alphaVector = linspace(0.6,alphaUpperLim,numOfAlpha);
alphaFB = 0.912; % 0.88 for 0.6-6, 0.95

channelDataStructFB = ChannelData(specFB,irfTFB,dtFB,1.5,1:size(specFB,2),[],1800);
Laguerre_StructFB = LaguerreModel(channelDataStructFB,LagOrder,alphaFB);
Laguerre_StructFB.estimate_laguerre([],0);
LTFB = Laguerre_StructFB.LTs;
plot(trueLTFB-LTFB)

fprintf('Finished FB\n')

%% convolusion V4
specV4 = filter(irfTV4,1,decayV4);
specV4 = specV4./max(specV4);
SNR = 55; %in dB, SNR = 20log10(Max/noise)
if SNR==0
    noise = zeros(size(specV4));
else
    noise = randn(size(specV4))*1/db2mag(SNR);
end
DC = 1:size(specV4,1);
DC = DC';
[~,MaxIdx] = max(specV4);
MaxIdx = mode(MaxIdx);
DC(DC<MaxIdx)=0;
DC(DC>=MaxIdx)=1;
DC = DC*0.0000;
specV4 = specV4+noise+DC;
figure
plot(irfTV4)
hold on
plot(specV4(:,1:50:end))
hold off
title('V4 Simulation data and irf')

%% deconvolution
alphaUpperLim=alpha_up(size(specV4,1),LagOrder,[],[])
% alphaUpperLim=0.916;

% alphaVector = linspace(0.6,alphaUpperLim,numOfAlpha);
alphaV4 = 0.9; % 0.88 for 0.6-6, 0.95

LaguerreFB = Laguerre_StructFB.LaguerreBasis;
LaguerreFB'*LaguerreFB
LaguerreV4 = interp1(tFB,LaguerreFB,tV4,'pchip');
T = LaguerreV4'*LaguerreV4
LaguerreV4 = LaguerreV4/sqrt(trace(T)/12);
LaguerreV4'*LaguerreV4
figure;plot(LaguerreV4)

K=11;
figure
plot(tFB,LaguerreFB(:,K),'r.')
hold on
plot(tV4,LaguerreV4(:,K),'go')


channelDataStructV4 = ChannelData(specV4,irfTV4,dtV4,1.5,1:size(specV4,2),[],1800);
Laguerre_StructV4 = LaguerreModel(channelDataStructV4,LagOrder,alphaV4);
figure;plot(Laguerre_StructV4.LaguerreBasis)
Laguerre_StructV4.LaguerreBasis = LaguerreV4;
Laguerre_StructV4.estimate_laguerre([],0);
LTV4 = Laguerre_StructV4.LTs;
fit = get(Laguerre_StructV4,'fit');

%%
WF = 2;
figure;plot(specV4(:,WF))
hold on
plot(fit(:,WF))
hold off
%%
figure
plot(trueLTV4-LTV4)
fprintf('Finished FB\n')
%% plot specific alpha
% figure('Position',[200 200 800 450]);
% scatter(trueLT,LTArray(:,end),'b.')
% hold on
% plot([0 20],[0 20],'r-')
% grid on
% xlim([0 20])
% ylim([0 20])
% xlabel('Ture lifetime (ns)')
% ylabel('Computed lifetime (ns)')
% title(sprintf('Tuncation = %.2f, Alpha = %.4f, SNR = %d, Laguerre order = %d',tWindow,alphaVector(end),SNR, LagOrder))
%% plot specific alpha
% figure('Position',[200 200 800 450]);
% scatter(trueLT,LTArray(:,end)-trueLT,'b.')
% xline(0.5,'--r','0.5 ns','LineWidth',1.5)
% yline(0,'--r','LineWidth',1.5)
% % hold on
% % plot([0 10],[0 10],'r-')
% grid on
% xlim([0 17])
% ylim([-2 2])
% xlabel('Ture lifetime (ns)')
% ylabel('Lifetime difference (ns)')
% title(sprintf('Tuncation = %.2f, Alpha = %.4f, SNR = %d, Laguerre order = %d',tWindow,alphaVector(end),SNR,LagOrder))
%% plot specific alpha
% figure('Position',[200 200 800 450]);
% scatter(tau1,LTArray(:,end)-trueLT,'b.')
% xline(0.5,'--r','0.5 ns','LineWidth',1.5)
% yline(0,'--r','LineWidth',1.5)
% % hold on
% % plot([0 10],[0 10],'r-')
% grid on
% xlim([0 17])
% ylim([-2 2])
% xlabel('tau1 lifetime (ns)')
% ylabel('Lifetime difference (ns)')
% title(sprintf('Tuncation = %.2f, Alpha = %.4f, SNR = %d, Laguerre order = %d',tWindow,alphaVector(end),SNR,LagOrder))
