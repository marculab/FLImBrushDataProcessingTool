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
dataObj = Ch2DataObj; %% set your channel here
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
LTLag = dataObj.LTsAll;
irf = dataObj.APDObj.irfTNorm;

%% exponential fit 
tic

order = 3; % set your order here
numOfDeconWF =  length(deconIdx);
ATemp = zeros(numOfDeconWF,order);
TauTemp = zeros(numOfDeconWF,order);
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
    [A, T, avglife, intensity, fit1Temp(:,idxSpecTemp), raw, decay1Temp(:,idxSpecTemp)] = multiexp_fit(specTempCh1,dataObj.dtUp,laser,order,[]);
    ATemp(idxSpecTemp,:) = A';
    TauTemp(idxSpecTemp,:) = T';
    LTexpTemp(idxSpecTemp) = avglife;
end
WTemp = ATemp./sum(ATemp,2);
res1Temp = WFAligned-fit1Temp;
% WTemp = WTemp';
%
A1 = zeros(numOfWF,order);
Tau1 = zeros(numOfWF,order);
W1 = zeros(numOfWF,order);
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

LTexp1Decay = h_lifet(decay1,dataObj.dtUp);

ExpResult1.A = A1;
ExpResult1.Tau = Tau1;
ExpResult1.LTexp = LTexp1;
ExpResult1.W = W1;
ExpResult1.Order = order;
ExpResult1.Fit = fit1;
ExpResult1.Decay = decay1;
ExpResult1.res = res1;

tausAvg = mean(Tau1)

processingT = toc;
fprintf('Exponential fitting time %.4f\n',processingT);
%% exponential fit with fixed tau
tic
ATemp = zeros(numOfDeconWF,order);
TauTemp = zeros(numOfDeconWF,order);
fit2Temp = zeros(size(WFAligned));
decay2Temp = zeros(size(WFAligned));

for i = 1:length(irfIdxValueCh1)
    idxTemp = irfIdxValueCh1(i);
    laser = irf(:,idxTemp);
    idxSpecTemp = find(irfIdx1==idxTemp);
    specTempCh1 = WFAligned(:,idxSpecTemp);
    [A, T, avglife, intensity, fit2Temp(:,idxSpecTemp), raw, decay2Temp(:,idxSpecTemp)] = multiexp_fit(specTempCh1,dataObj.dtUp,laser,order,tausAvg);
    ATemp(idxSpecTemp,:) = A';
    TauTemp(idxSpecTemp,:) = T';
    LTexpTemp(idxSpecTemp) = avglife;
end
WTemp = ATemp./sum(ATemp,2);
res2Temp = WFAligned-fit2Temp;
% WTemp = WTemp';
%
A2 = zeros(numOfWF,order);
Tau2 = zeros(numOfWF,order);
W2 = zeros(numOfWF,order);
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

LTexp2Decay = h_lifet(decay2,dataObj.dtUp);


ExpResult2.A = A2;
ExpResult2.Tau = Tau2;
ExpResult2.LTexp = LTexp2;
ExpResult2.W = W2;
ExpResult2.Order = order;
ExpResult2.Fit = fit2;
ExpResult2.Decay = decay2;
ExpResult2.res = res2;

processingT = toc;
fprintf('Exponential fitting with fixed taus time %.4f\n',processingT);



%% plot Laguerre
idx = round(rand*numOfWF);

figure('Position',[300 300 1200 400])
tiledlayout(1,3)
nexttile
plot(WFALignedAll(:,idx),'b.','MarkerSize',10)
hold on
plot(LagFittingAll(:,idx),'r-','LineWidth',1)
plot(LagResAll(:,idx)-0.1,'m.')
yline(-0.1+0.025,'r--','res=0.025')
yline(-0.1-0.025,'r--')
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
yticks((-0.2:0.1:1))
yticklabels({'','res=0','0','0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','1.0',})
grid on
grid minor
axis tight
ylim([-0.2 1])
legend('raw data', 'fitting', 'residue')
title(sprintf('Exponential fitting, order = %d, WF index %d,\n lifetime = %.3f(formula), lifetime = %.3f(decay) \n a1=%.3f, a2=%.3f, a3=%.3f, tau1=%.3f, tau2=%.3f, tau3=%.3f', ...
    order,idx, LTexp1(idx), LTexp1Decay(idx), W1(idx,1), W1(idx,2), W1(idx,3), Tau1(idx,1), Tau1(idx,2), Tau1(idx,3)))


% plot tri-exponential
nexttile
plot(WFALignedAll(:,idx),'b.','MarkerSize',10)
hold on
plot(fit2(:,idx),'r-','LineWidth',1)
plot(res2(:,idx)-0.1,'m.')
yline(-0.1+0.025,'r--','res=0.025')
yline(-0.1-0.025,'r--')
yticks((-0.2:0.1:1))
yticklabels({'','res=0','0','0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','1.0',})
grid on
grid minor
axis tight
ylim([-0.2 1])
legend('raw data', 'fitting', 'residue')
title(sprintf('Exponential fitting (fixed taus), order = %d, WF index %d,\n lifetime = %.3f(formula), lifetime = %.3f(decay) \n a1=%.3f, a2=%.3f, a3=%.3f, tau1=%.3f, tau2=%.3f, tau3=%.3f', ...
    order,idx, LTexp2(idx), LTexp2Decay(idx), W2(idx,1), W2(idx,2), W2(idx,3), Tau2(idx,1), Tau2(idx,2),Tau2(idx,3)))

