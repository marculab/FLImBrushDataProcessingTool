% code to run simulated deconvolution to determine the best alpha value

%%
clear all
close all
clc

%% get random lifetime
N = 2500;
% trueLT = rand(N,1);
lowLT = 0.3;
highLT = 15;
trueLT = lowLT+rand(N,1)*(highLT-lowLT);
binEdge = linspace(lowLT,highLT,40);
figure
histogram(trueLT,binEdge);
title('histogram of randomly distributed lifetimes')
xlabel('Lifetime (ns)')
ylabel('Count')

%% get decays
tWindow = 80; % ns same for V4 nad V5
dt = 0.2;
t = 0:dt:tWindow-dt;
decay = zeros(tWindow/dt,N);
for i = 1:N
    tau = trueLT(i);
    decay(:,i) = exp(-t/tau);
end

plot(decay(:,1:50:end))

%% load and truncate irf
load('Ch1iRF.mat')
[~,irfMaxIdx] = max(irf);
tLength = tWindow/dt;
startIdx = round(irfMaxIdx-0.1*tLength);
if startIdx <=0
    startIdx = 1;
end
irfT = irf(startIdx:startIdx+tLength-1);
figure
plot(irfT)
title('irf')
%%
spec = filter(irfT,1,decay);
figure
plot(irfT)
hold on
plot(spec(:,1:50:end))
hold off
title('Simulation data and irf')

%% deconvolution
alphaUpperLim=alpha_up(size(spec,1),12,[],[]);
numOfAlpha = 50;
alphaVector = linspace(0.5,alphaUpperLim,numOfAlpha);
LTArray = zeros(N,numOfAlpha);
f = waitbar(0,'Starting');
LagOrder = 20;
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
ylim([0 15])
zlim([0 15])
xlabel('Alpha')
ylabel('True Lifetime (ns)')
zlabel('Computed Lifetime (ns)')
title(sprintf('Order = 12, time window = %dns, %d simulation per alpha',tWindow,N));
saveas(gcf,'Lifetime Plot.fig')
%% plot 3D error
E = trueLT-LTArray;
figure
scatter3(XX(:),YY(:),E(:),'b.');
ylim([0 15])
% zlim([0 15])
xlabel('Alpha')
ylabel('True Lifetime (ns)')
zlabel('Error (ns)')
title(sprintf('Order = %d, time window = %dns, %d simulation per alpha',LagOrder,tWindow,N));
saveas(gcf,'Lifetime Absolute Error Plot.fig')

%% plot 3D % error
E = trueLT-LTArray;
EE = round(E./trueLT*10000)/100;
figure
scatter3(XX(:),YY(:),EE(:),'b.');
% surf(XX(:),YY(:),EE(:));
ylim([0 15])
zlim([-20 50])
xlabel('Alpha')
ylabel('True Lifetime (ns)')
zlabel('Percentage Error')
title(sprintf('Order = 12, time window = %dns, %d simulation per alpha',tWindow,N));
saveas(gcf,'Lifetime Relative Error Plot.fig')
