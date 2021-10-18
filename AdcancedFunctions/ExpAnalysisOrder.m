% code to run multiexpoential fit to FLImBRUSH output

%%
clear
close all
clc

%% add path
% addpath(genpath('./MultiExp Generic')) %add path to exponential code
addpath(genpath('../Algorithms'))

%% select file
root = 'D:\LoaclData\5ALAFLImTest\Subject032_20210811';
[DeConfile,DeConpath] = uigetfile([root '\*.mat'],'Please select DeCon file','MultiSelect','off');

%% load in file
P = fullfile(DeConpath,DeConfile);
load(P)

%% set up channel
dataObj = Ch2DataObj; %% set up your channel here
%% get Laguerre result
WFAligned = get(dataObj,'wf_aligned');
LagFitting = get(dataObj,'fit');
LagRes = get(dataObj,'res');
numOfWF = dataObj.numOfAvgWFs; % get number of waveforms
deconIdx = dataObj.deconIdx;
WFALignedAll = zeros(size(WFAligned,1),numOfWF);
LagFittingAll = zeros(size(LagFitting,1),numOfWF);
LagResAll = zeros(size(LagRes,1),numOfWF);
WFALignedAll(:,deconIdx) = WFAligned;
LagFittingAll(:,deconIdx) = LagFitting;
LagResAll(:,deconIdx) = LagRes;
LagSE = sum(LagResAll.^2);
LTLag = dataObj.LTsAll;
irf = dataObj.APDObj.irfTNorm;

%% bi-exponential fit
order1 = 2;
numOfDeconWF =  length(deconIdx);
ATemp = zeros(numOfDeconWF,order1);
TauTemp = zeros(numOfDeconWF,order1);
LTexpTemp = zeros(numOfDeconWF,1);
irfIdx1 = dataObj.irfIdx;
irfIdxValueCh1 = unique(irfIdx1,'sorted');
fit1Temp = zeros(size(WFAligned));
decay1Temp = zeros(size(WFAligned));

for i = 1:length(irfIdxValueCh1)
    idxTemp = irfIdxValueCh1(i);
    laser = irf(:,idxTemp);
    idxSpecTemp = find(irfIdx1==idxTemp);
    specTempCh1 = WFAligned(:,idxSpecTemp);
    [A, T, avglife, intensity, fit1Temp(:,idxSpecTemp), raw, decay1Temp(:,idxSpecTemp)] = multiexp_fit(specTempCh1,dataObj.dtUp,laser,order1,[]);
    ATemp(idxSpecTemp,:) = A';
    TauTemp(idxSpecTemp,:) = T';
    LTexpTemp(idxSpecTemp) = avglife;
end
WTemp = ATemp./sum(ATemp,2);
res1Temp = WFAligned-fit1Temp;
% WTemp = WTemp';
%
A1 = zeros(numOfWF,order1);
Tau1 = zeros(numOfWF,order1);
W1 = zeros(numOfWF,order1);
LTexp1 = zeros(numOfWF,1);
fit1 = zeros(size(fit1Temp,1),numOfWF);
decay1 = zeros(size(decay1Temp,1),numOfWF);
res1 = zeros(size(decay1Temp,1),numOfWF);

A1(deconIdx,:) = ATemp;
Tau1(deconIdx,:) = TauTemp;
W1(deconIdx,:) = WTemp;
LTexp1(deconIdx) = LTexpTemp;
fit1(:,deconIdx) = fit1Temp;
decay1(:,deconIdx) = decay1Temp;
res1(:,deconIdx) = res1Temp;

res1SE = sum(res1.^2);
LTexp1Decay = h_lifet(decay1,dataObj.dtUp);

ExpResult1.A = A1;
ExpResult1.Tau = Tau1;
ExpResult1.LTexp = LTexp1;
ExpResult1.W = W1;
ExpResult1.Order = order1;
ExpResult1.Fit = fit1;
ExpResult1.Decay = decay1;
ExpResult1.res = res1;

taus1Avg = mean(Tau1)

%% tri-exponential fit
order2 = 3;
ATemp = zeros(numOfDeconWF,order2);
TauTemp = zeros(numOfDeconWF,order2);
fit2Temp = zeros(size(WFAligned));
decay2Temp = zeros(size(WFAligned));

for i = 1:length(irfIdxValueCh1)
    idxTemp = irfIdxValueCh1(i);
    laser = irf(:,idxTemp);
    idxSpecTemp = find(irfIdx1==idxTemp);
    specTempCh1 = WFAligned(:,idxSpecTemp);
    [A, T, avglife, intensity, fit2Temp(:,idxSpecTemp), raw, decay2Temp(:,idxSpecTemp)] = multiexp_fit(specTempCh1,dataObj.dtUp,laser,order2,[]);
    ATemp(idxSpecTemp,:) = A';
    TauTemp(idxSpecTemp,:) = T';
    LTexpTemp(idxSpecTemp) = avglife;
end
WTemp = ATemp./sum(ATemp,2);
res2Temp = WFAligned-fit2Temp;
% WTemp = WTemp';
%
A2 = zeros(numOfWF,order2);
Tau2 = zeros(numOfWF,order2);
W2 = zeros(numOfWF,order2);
LTexp2 = zeros(numOfWF,1);
fit2 = zeros(size(fit2Temp,1),numOfWF);
decay2 = zeros(size(decay2Temp,1),numOfWF);
res2 = zeros(size(fit2Temp,1),numOfWF);


A2(deconIdx,:) = ATemp;
Tau2(deconIdx,:) = TauTemp;
W2(deconIdx,:) = WTemp;
LTexp2(deconIdx) = LTexpTemp;
fit2(:,deconIdx) = fit2Temp;
decay2(:,deconIdx) = decay2Temp;
res2(:,deconIdx) = res2Temp;

res2SE = sum(res2.^2);
LTexp2Decay = h_lifet(decay2,dataObj.dtUp);


ExpResult2.A = A2;
ExpResult2.Tau = Tau2;
ExpResult2.LTexp = LTexp2;
ExpResult2.W = W2;
ExpResult2.Order = order2;
ExpResult2.Fit = fit2;
ExpResult2.Decay = decay2;
ExpResult2.res = res2;

taus2Avg = mean(Tau2)

%% tri-exponential fit
order4 = 4;
ATemp = zeros(numOfDeconWF,order4);
TauTemp = zeros(numOfDeconWF,order4);
fit4Temp = zeros(size(WFAligned));
decay4Temp = zeros(size(WFAligned));

for i = 1:length(irfIdxValueCh1)
    idxTemp = irfIdxValueCh1(i);
    laser = irf(:,idxTemp);
    idxSpecTemp = find(irfIdx1==idxTemp);
    specTempCh1 = WFAligned(:,idxSpecTemp);
    [A, T, avglife, intensity, fit4Temp(:,idxSpecTemp), raw, decay4Temp(:,idxSpecTemp)] = multiexp_fit(specTempCh1,dataObj.dtUp,laser,order4,[]);
    ATemp(idxSpecTemp,:) = A';
    TauTemp(idxSpecTemp,:) = T';
    LTexpTemp(idxSpecTemp) = avglife;
end
WTemp = ATemp./sum(ATemp,2);
res4Temp = WFAligned-fit4Temp;
% WTemp = WTemp';
%
A4 = zeros(numOfWF,order4);
Tau4 = zeros(numOfWF,order4);
W4 = zeros(numOfWF,order4);
LTexp4 = zeros(numOfWF,1);
fit4 = zeros(size(fit4Temp,1),numOfWF);
decay4 = zeros(size(decay4Temp,1),numOfWF);
res4 = zeros(size(fit4Temp,1),numOfWF);


A4(deconIdx,:) = ATemp;
Tau4(deconIdx,:) = TauTemp;
W4(deconIdx,:) = WTemp;
LTexp4(deconIdx) = LTexpTemp;
fit4(:,deconIdx) = fit4Temp;
decay4(:,deconIdx) = decay4Temp;
res4(:,deconIdx) = res4Temp;

res4SE = sum(res4.^2);
LTexp4Decay = h_lifet(decay4,dataObj.dtUp);

ExpResult4.A = A4;
ExpResult4.Tau = Tau4;
ExpResult4.LTexp = LTexp4;
ExpResult4.W = W4;
ExpResult4.Order = order4;
ExpResult4.Fit = fit4;
ExpResult4.Decay = decay4;
ExpResult4.res = res4;

taus4Avg = mean(Tau4)
%% plot Laguerre
idx = round(rand*numOfWF);

figure('Position',[300 300 1200 400])
tiledlayout(2,2)
nexttile
plot(WFALignedAll(:,idx),'b.','MarkerSize',10)
hold on
plot(LagFittingAll(:,idx),'r-','LineWidth',1)
plot(LagResAll(:,idx)-0.1,'m.')
yline(-0.1+0.025,'r--','res=0.025')
yline(-0.1-0.025,'r--')
text(150, -0.15, ['\Sigma res = ', num2str(LagSE(idx))])
yticks((-0.2:0.1:1))
yticklabels({'','res=0','0','0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','1.0',})
grid on
grid minor
axis tight
ylim([-0.2 1])
legend('raw data', 'fitting', 'residue')
title(sprintf('Laguerre fitting, WF index %d,\n lifetime = %.3f',idx, LTLag(idx)))

% plot bi-exponential
nexttile
plot(WFALignedAll(:,idx),'b.','MarkerSize',10)
hold on
plot(fit1(:,idx),'r-','LineWidth',1)
plot(res1(:,idx)-0.1,'m.')
yline(-0.1+0.025,'r--','res=0.025')
yline(-0.1-0.025,'r--')
text(150, -0.15, ['\Sigma res = ', num2str(res1SE(idx))])
yticks((-0.2:0.1:1))
yticklabels({'','res=0','0','0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','1.0',})
grid on
grid minor
axis tight
ylim([-0.2 1])
legend('raw data', 'fitting', 'residue')
title(sprintf('Exponential fitting, order = %d, WF index %d,\n lifetime = %.3f(formula), lifetime = %.3f(decay),\n a1=%.3f, a2=%.3f, tau1=%.3f, tau2=%.3f', ...
    order1,idx, LTexp1(idx), LTexp1Decay(idx), W1(idx,1), W1(idx,2), Tau1(idx,1), Tau1(idx,2)))

% plot tri-exponential
nexttile
plot(WFALignedAll(:,idx),'b.','MarkerSize',10)
hold on
plot(fit2(:,idx),'r-','LineWidth',1)
plot(res2(:,idx)-0.1,'m.')
yline(-0.1+0.025,'r--','res=0.025')
yline(-0.1-0.025,'r--')
text(150, -0.15, ['\Sigma res = ', num2str(res2SE(idx))])
yticks((-0.2:0.1:1))
yticklabels({'','res=0','0','0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','1.0',})
grid on
grid minor
axis tight
ylim([-0.2 1])
legend('raw data', 'fitting', 'residue')
title(sprintf('Exponential fitting, order = %d, WF index %d,\n lifetime = %.3f(formula), lifetime = %.3f(decay) \n a1=%.3f, a2=%.3f, a3=%.3f, tau1=%.3f, tau2=%.3f, tau3=%.3f', ...
    order2,idx, LTexp2(idx), LTexp2Decay(idx), W2(idx,1), W2(idx,2), W2(idx,3), Tau2(idx,1), Tau2(idx,2),Tau2(idx,3)))

% plot quad-exponential
nexttile
plot(WFALignedAll(:,idx),'b.','MarkerSize',10)
hold on
plot(fit4(:,idx),'r-','LineWidth',1)
plot(res4(:,idx)-0.1,'m.')
yline(-0.1+0.025,'r--','res=0.025')
yline(-0.1-0.025,'r--')
text(150, -0.15, ['\Sigma res = ', num2str(res4SE(idx))])
yticks((-0.2:0.1:1))
yticklabels({'','res=0','0','0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','1.0',})
grid on
grid minor
axis tight
ylim([-0.2 1])
legend('raw data', 'fitting', 'residue')
title(sprintf('Exponential fitting, order = %d, WF index %d,\n lifetime = %.3f(formula), lifetime = %.3f(decay) \n a1=%.3f, a2=%.3f, a3=%.3f, a4=%.3f, tau1=%.3f, tau2=%.3f, tau3=%.3f, tau4=%.3f', ...
    order4,idx, LTexp4(idx), LTexp4Decay(idx), W4(idx,1), W4(idx,2), W4(idx,3), W4(idx,4), Tau4(idx,1), Tau4(idx,2), Tau4(idx,3), Tau4(idx,4)))

%% plot square error
figure
plot(LagSE,'.','Color','#0072BD')
hold on
plot(res1SE,'.','Color','#D95319')
plot(res2SE,'.','Color','#77AC30')
plot(res4SE,'.','Color','#A2142F')
set(gca, 'YScale', 'log')
axis tight
grid on
grid minor
% ylim([0 0.05])
legend('Laguerre','exp 2','exp 3','exp 4')
xlabel('Waveform index')
ylabel('Square Error')
