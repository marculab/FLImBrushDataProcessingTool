% code to batch reprocess FLImBRUSH .mat data
% Xiangnan 09/09/2022

%% clear workspace
clear
close all
clc

%% load processed .mat file
load('C:\Users\Xiangnan\Box Sync\FLImBrush vs V4\20210722V4vsFLImBRUSHDye\FLImBRUSH\20210722DyeTest2\Decon 2022-8-30 70ns 0.916\20210722DyeTest2_03.mat')
%% copy old stucture
Ch1DataObj_old = Ch1DataObj;
Ch2DataObj_old = Ch2DataObj;
Ch3DataObj_old = Ch3DataObj;
EOP_H1G_old = EOP_H1G;
EOP_H1S_old = EOP_H1S;
SP_G_old = SP_G;
SP_S_old = SP_S;

clear('Ch1DataObj','Ch2DataObj','Ch3DataObj','EOP_H1G','EOP_H1S','SP_G','SP_S')
%% reload BG to new sampling rate
upSampleFactor = 5;
bgObj = backGround(fullfile(dataInfoObj.bgFileFolder,dataInfoObj.bgFileName));
loadBG(bgObj,upSampleFactor);
%% retrive common data to all channels
laserRepRateIn = Ch1DataObj_old.laserRepRate;
LVAvgIn = Ch1DataObj_old.dataAveragingLV;
bgLow = Ch1DataObj_old.bgLow/4*5;
bgLow = round(bgLow);
bgHigh = Ch1DataObj_old.bgHigh/4*5;
bgHigh = round(bgHigh);
avg = Ch1DataObj_old.dataAveraging;
dataLow = Ch1DataObj_old.dataLow;
dataHigh = Ch1DataObj_old.dataHigh;
Truncation = 680*0.08; % same as V4
exclude = Ch1DataObj_old.exclude;
LaguerreOrder = 12;
alpha = 0.916;
expOrder = Ch1DataObj_old.expDeconObj.order;
%% retrive nessesery data for reprocessing
Ch1rawData = Ch1DataObj_old.rawData;
Ch1CtrlV = Ch1DataObj_old.rawCtrlV;
Ch1APDObjIn = Ch1DataObj_old.APDObj;
Ch1dtIn = Ch1DataObj_old.dtRaw;

%% update channel data obj to new decon setting
figure; tiledlayout('flow') ;
Ch1DataObj = ChannelDataAPD(Ch1rawData, Ch1CtrlV, LVAvgIn, Ch1APDObjIn, Ch1dtIn, bgObj.bgCh1, laserRepRateIn, upSampleFactor);
nexttile
plot(Ch1rawData(:,1:500:end))
title('Raw Data')

removeDCData(Ch1DataObj);
nexttile
plot(Ch1DataObj.rawDataDCRemoved(:,1:500:end));
title('DC removed Raw Data')

upSampleData(Ch1DataObj);
nexttile
plot(Ch1DataObj.rawDataUpsampled(:,1:500:end));
title('Upsampled Data')

alignWF_CFD(Ch1DataObj, 1, 1:680)
removeSaturation(Ch1DataObj,dataLow,dataHigh);
avgData(Ch1DataObj, avg)
nexttile
plot(Ch1DataObj.averagedData(:,1:500:end));
title('Averaged Data')

removeDCBG(Ch1DataObj,bgLow,bgHigh,0.03);
nexttile
plot(Ch1DataObj.preProcessedData(:,1:500:end));
title('BG removed Data')

truncateData(Ch1DataObj,Truncation);
nexttile
plot(Ch1DataObj.dataT(:,1:500:end));
title('Truncated Data')

runDeconLG(Ch1DataObj,exclude,LaguerreOrder,alpha);

runDeconExp(Ch1DataObj,expOrder,[],[],[],exclude);

%% plot fitting
Idx = round(rand()*Ch1DataObj.numOfAvgWFs);
rawDataToPlot = get(Ch1DataObj,'wf_aligned',Idx);
rawDataMax = max(rawDataToPlot);
fitToPlotLG = get(Ch1DataObj,'fit',Idx);
fitToPlotExp = get(Ch1DataObj.expDeconObj,'fit',Idx);
residueToPlotLG = rawDataToPlot-fitToPlotLG;
residueToPlotExp = rawDataToPlot-fitToPlotExp;
gain = Ch1DataObj.gain(Idx);
CtrlV = Ch1DataObj.CtrlV(Idx);
LT_LG =  Ch1DataObj.Lg_LTs(Idx);
irfIdx = Ch1DataObj.irfIdx(Idx);
irf = Ch1DataObj.APDObj.irfTNorm(:,irfIdx);
t = Ch1DataObj.dtUp*(1:length(rawDataToPlot))';
rawColor = [0 0 1 0.5];

figure
plot(t,rawDataToPlot,'.-','Color', rawColor,'LineWidth', 1, 'MarkerSize',10)
hold on
plot(t(1:length(irf)),rawDataMax/max(irf)*irf,'Color',[0 0 0 0.5],'LineWidth', 1)
plot(t,fitToPlotLG,'m','LineWidth', 1.2)
plot(t,fitToPlotExp,'g-.','LineWidth', 1.2)
hold off
lgd = legend('Raw Data', 'iRF', 'LG fit','mExp fit');
lgd.FontSize = 8;
title_temp = sprintf('Point %2d, Laguerre Lifetime %3.2f, Gain %4.0f.',Idx,LT_LG,gain);
title_temp = [1 title_temp];
title(title_temp)
% ylim([-0.2 Ch1DataObj.dataHigh+0.1])

%%
LifetimeOld = Ch1DataObj_old.Lg_LTs;
Lifetime = Ch1DataObj.Lg_LTs;
figure
tiledlayout(1,2)
nexttile
scatter(Lifetime,Lifetime-LifetimeOld)
nexttile
histogram(Lifetime-LifetimeOld)