% code to run simulated deconvolution to determine the best alpha value

%% clear workspace and add code to MATLAB path
clear 
close all
clc
addpath(genpath('C:\Users\Xiangnan\Documents\MyGitRepo\FLImBrushDataProcessingTool'))
%% get random lifetime
N = 10000; %# of random wavefoms
% trueLT = rand(N,1);
tau1Low = 0.3;
tau1High = 4;
tau2 = 6;
tau1 = tau1Low+rand(N,1)*(tau1High-tau1Low);

k=9; % photon ratio

%% get decays
dt = 0.1; % time resolution
tWindow = 75; % total data length in time (ns)
trueLTC = zeros(size(tau1));
t = 0:dt:tWindow-dt;
decay = zeros(round(tWindow/dt),N); % decays
for i = 1:N % loop through all lifetimes
    M = [1 1;tau1(i) -k*tau2];
    A = M\[1;0]; % get pre-exponential factor
    trueLTC(i) = (A(1)*tau1(i).^2+A(2)*tau2^2)./(A(1)*tau1(i)+A(2)*tau2); % true lifetime from formula
    decay(:,i) = A(1)*exp(-t/tau1(i))+A(2)*exp(-t/tau2);
    
end
trueLT = h_lifet(decay,dt,'average')'; % true lifetime from decay
binEdge = linspace(min(trueLT),max(trueLT),40); % bin edge for histogram
%------------------------------plot-----------------------------------
figure
scatter(tau1,trueLTC-trueLT);
% title('histogram of randomly distributed lifetimes')
xlabel('tau1 (ns)')
ylabel('Lifetime difference')

figure
scatter(tau1,trueLTC);
title('Average lifetime vs tau1, fixed tau2 = 16 ns')
xlabel('tau1 (ns)')
ylabel('Lifetime')

figure
plot(decay(:,1:500:end))

%% load and truncate V4 iIRF
irfStruc = load('TriplexV4_IRF_H&N_Data_Analysis.mat'); % load iIRF
irf = irfStruc.laser{1}; % get channel 1 iIRF
irfT = irf;
figure
plot(irfT)
title('irf')
%% simulate waveforms
spec = filter(irfT,1,decay); % convolution
spec = spec./max(spec); % normalize waveform
SNR = 0; % SNR in dB, SNR = 20log10(Max/noise)
if SNR==0
    noise = zeros(size(spec));
else
    noise = randn(size(spec))*1/db2mag(SNR);
end
DC = 1:size(spec,1); % DC value, this is the verticle shift due to digitizer
DC = DC';
% [~,MaxIdx] = max(spec);
% MaxIdx = mode(MaxIdx);
% DC(DC<MaxIdx)=0;
% DC(DC>=MaxIdx)=1;
DC = DC*0.0000; % update DC value
spec = spec+noise+DC;
figure
plot(irfT)
hold on
plot(spec(:,1:50:end))
hold off
title('Simulation data and irf')

%% deconvolution, you can replace this by the version of code that does not use the objects
LagOrder = 12;
alphaUpperLim=alpha_up(size(spec,1),LagOrder,[],[]);
% alphaUpperLim=0.916;

numOfAlpha = 1;
% alphaVector = linspace(0.6,alphaUpperLim,numOfAlpha);
alphaVector = 0.910; % 0.88 for 0.6-6, 0.95
LTArray = zeros(N,numOfAlpha);
f = waitbar(0,'Starting');
for i=1:numOfAlpha
    alphaTemp = alphaVector(i);
    channelDataStruct = ChannelData(spec,irfT,dt,1.5,1:size(spec,2),[],1800);
    Laguerre_Struct = LaguerreModel(channelDataStruct,LagOrder,alphaTemp);
    Laguerre_Struct.estimate_laguerre(0,0);
    LTArray(:,i) = Laguerre_Struct.LTs;
    waitbar(i/numOfAlpha,f,sprintf('Progress: %d %%',round(i/numOfAlpha*100)));
end
fprintf('Finished\n')
delete(f)
%% repopulate array for plotting
XX = repmat(alphaVector,[N,1]); % alpha
YY = repmat(trueLT,[1,numOfAlpha]); %$ true lifetime
%% plot 3D 
figure
scatter3(XX(:),YY(:),LTArray(:),'b.');
ylim([0 20])
zlim([0 20])
xlabel('Alpha')
ylabel('True Lifetime (ns)')
zlabel('Computed Lifetime (ns)')
title(sprintf('Order = %d, time window = %f ns, %d simulation per alpha',LagOrder,tWindow,N));
% saveas(gcf,'Lifetime Plot.fig')
%% plot 3D error
% E = LTArray-trueLT;
% figure
% scatter3(XX(:),YY(:),E(:),'b.');
% ylim([0 10])
% zlim([-2 2])
% xlabel('Alpha')
% ylabel('True Lifetime (ns)')
% zlabel('Error (ns)')
% title(sprintf('Order = %d, time window = %f ns, %d simulation per alpha',LagOrder,tWindow,N));
% saveas(gcf,'Lifetime Absolute Error Plot.fig')

%% plot 3D % error
% E = trueLT-LTArray;
% EE = round(E./trueLT*10000)/100;
% figure
% scatter3(XX(:),YY(:),EE(:),'b.');
% % surf(XX(:),YY(:),EE(:));
% ylim([0 10])
% zlim([-2 2])
% xlabel('Alpha')
% ylabel('True Lifetime (ns)')
% zlabel('Percentage Error')
% title(sprintf('Order = %d, time window = %f ns, %d simulation per alpha',LagOrder,tWindow,N));
% saveas(gcf,'Lifetime Relative Error Plot.fig')

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
