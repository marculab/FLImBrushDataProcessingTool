% code to run simulated deconvolution to determine the best alpha value

%%
clear 
close all
clc

%% get random lifetime
N = 10000;
% trueLT = rand(N,1);
lowLT = 0.3;
highLT = 20;
trueLT = lowLT+rand(N,1)*(highLT-lowLT);
binEdge = linspace(lowLT,highLT,40);
figure
histogram(trueLT,binEdge);
title('histogram of randomly distributed lifetimes')
xlabel('Lifetime (ns)')
ylabel('Count')

%% get decays
dt = 0.1;
tWindow = 154; % ns same for V4 nad V5

t = 0:dt:tWindow-dt;
decay = zeros(round(tWindow/dt),N);
for i = 1:N
    tau = trueLT(i);
    decay(:,i) = exp(-t/tau);
end

% plot(decay(:,1:50:end))
figure
plot(t,decay(:,5))
grid on
xlabel('Time (ns)')
ylabel('Amplitude')
title(['Lifetime ' num2str(trueLT(5))])

%% load and truncate irf
irfStruc = load('..\APDDetectorFile\M00549707_DCS.mat');
UpFactor = irfStruc.irfRawdt/dt;
irf = interp(irfStruc.irf(:,120),UpFactor);
[~,irfMaxIdx] = max(irf);
tLength = tWindow/dt;
startIdx = round(irfMaxIdx-0.1*tLength);
if startIdx <=0
    startIdx = 1;
end
try
irfT = irf(startIdx:startIdx+tLength-1);
catch
    irfT = irf(startIdx:end);
end
figure
plot(irfT)
title('irf')
%%
spec = filter(irfT,1,decay);
spec = spec./max(spec);
ref = spec*0.000;
ref = circshift(ref, 32/dt);
SNR = 50; %in dB, SNR = 20log10(Max/noise)
if SNR==0
    noise = zeros(size(spec));
else
    noise = randn(size(spec))*1/db2mag(SNR);
end
DC = 1:size(spec,1);
DC = DC';
[~,MaxIdx] = max(spec);
MaxIdx = mode(MaxIdx);
DC(DC<MaxIdx)=0;
DC(DC>=MaxIdx)=1;
DC = DC*0.0000;
spec = spec+noise+DC+ref;
figure
plot(irfT)
hold on
plot(spec(:,1:50:end))
hold off
title('Simulation data and irf')

%% deconvolution
alphaUpperLim=alpha_up(size(spec,1),12,[],[]);
% alphaUpperLim=0.916;

numOfAlpha = 1;
% alphaVector = linspace(0.6,alphaUpperLim,numOfAlpha);
alphaVector = 0.8; % 0.88 for 0.6-6, 0.95
LTArray = zeros(N,numOfAlpha);
f = waitbar(0,'Starting');
LagOrder = 12;
for i=1:numOfAlpha
    alphaTemp = alphaVector(i);
    channelDataStruct = ChannelData(spec,irfT,dt,1.5,1:size(spec,2),[],1800);
    Laguerre_Struct = LaguerreModel(channelDataStruct,LagOrder,alphaTemp);
    Laguerre_Struct.estimate_laguerre(0);
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
figure('Position',[200 200 800 450]);
scatter(trueLT,LTArray(:,end),'b.')
hold on
plot([0 20],[0 20],'r-')
grid on
xlim([0 20])
ylim([0 20])
xlabel('Ture lifetime (ns)')
ylabel('Computed lifetime (ns)')
title(sprintf('Tuncation = %.2f, Alpha = %.4f, SNR = %d',tWindow,alphaVector(end),SNR))
%% plot specific alpha

figure('Position',[200 200 800 450]);
scatter(trueLT,LTArray(:,end)-trueLT,'b.')
yline(0,'--r','LineWidth',1.5)
% hold on
% plot([0 10],[0 10],'r-')
grid on
xlim([0 17])
ylim([-1 1])
xlabel('Ture lifetime (ns)')
ylabel('Lifetime difference (ns)')
title(sprintf('Tuncation = %.2f, Alpha = %.4f, SNR = %d',tWindow,alphaVector(end),SNR))

