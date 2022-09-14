% code to batch reprocess FLImBRUSH .mat data
% Xiangnan 09/09/2022
function ReprocessFB(input_mat_name)
% input_mat_name must be absolute path

[filepath,name,~] = fileparts(input_mat_name);
dataRoot = fileparts(filepath);
bgFolder = fullfile(dataRoot,'Background'); % find background folder as the data is moved
%% load processed .mat file
temp = load(input_mat_name,'dataInfoObj','Ch1DataObj','Ch2DataObj','Ch3DataObj');
dataInfoObj = temp.dataInfoObj;
Ch1DataObj_old = temp.Ch1DataObj;
Ch2DataObj_old = temp.Ch2DataObj;
Ch3DataObj_old = temp.Ch3DataObj;

%% reload BG to new sampling rate
upSampleFactor = 5;
bgObj = backGround(fullfile(bgFolder,dataInfoObj.bgFileName));
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
%% retrive nessesery Ch1 data for reprocessing
Ch1rawData = Ch1DataObj_old.rawData;
Ch1CtrlV = Ch1DataObj_old.rawCtrlV;
Ch1APDObjIn = Ch1DataObj_old.APDObj;
Ch1dtIn = Ch1DataObj_old.dtRaw;

% update channel1 data obj to new decon setting
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

alignWF_CFD(Ch1DataObj, 1, (50:170)*upSampleFactor)
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


%% retrive nessesery Ch2 data for reprocessing
plotStep = 100;
Ch2rawData = Ch2DataObj_old.rawData;
Ch2CtrlV = Ch2DataObj_old.rawCtrlV;
Ch2APDObjIn = Ch2DataObj_old.APDObj;
Ch2dtIn = Ch2DataObj_old.dtRaw;

% update channel1 data obj to new decon setting
figure; tiledlayout('flow') ;
Ch2DataObj = ChannelDataAPD(Ch2rawData, Ch2CtrlV, LVAvgIn, Ch2APDObjIn, Ch2dtIn, bgObj.bgCh2, laserRepRateIn, upSampleFactor);
nexttile
plot(Ch2rawData(:,1:plotStep:end))
title('Raw Data')

removeDCData(Ch2DataObj);
nexttile
plot(Ch2DataObj.rawDataDCRemoved(:,1:plotStep:end));
title('DC removed Raw Data')

upSampleData(Ch2DataObj);
nexttile
plot(Ch2DataObj.rawDataUpsampled(:,1:plotStep:end));
title('Upsampled Data')

alignWF_CFD(Ch2DataObj, 0.5)
removeSaturation(Ch2DataObj,dataLow,dataHigh);
avgData(Ch2DataObj, avg)
nexttile
plot(Ch2DataObj.averagedData(:,1:plotStep:end));
title('Averaged Data')

removeDCBG(Ch2DataObj,bgLow,bgHigh,0.03);
nexttile
plot(Ch2DataObj.preProcessedData(:,1:plotStep:end));
title('BG removed Data')

truncateData(Ch2DataObj,Truncation);
nexttile
plot(Ch2DataObj.dataT(:,1:plotStep:end));
title('Truncated Data')

runDeconLG(Ch2DataObj,exclude,LaguerreOrder,alpha);

runDeconExp(Ch2DataObj,expOrder,[],[],[],exclude);


%% retrive nessesery Ch3 data for reprocessing
plotStep = 100;
Ch3rawData = Ch3DataObj_old.rawData;
Ch3CtrlV = Ch3DataObj_old.rawCtrlV;
Ch3APDObjIn = Ch3DataObj_old.APDObj;
Ch3dtIn = Ch3DataObj_old.dtRaw;

% update channel1 data obj to new decon setting
figure; tiledlayout('flow') ;
Ch3DataObj = ChannelDataAPD(Ch3rawData, Ch3CtrlV, LVAvgIn, Ch3APDObjIn, Ch3dtIn, bgObj.bgCh3, laserRepRateIn, upSampleFactor);
nexttile
plot(Ch3rawData(:,1:plotStep:end))
title('Raw Data')

removeDCData(Ch3DataObj);
nexttile
plot(Ch3DataObj.rawDataDCRemoved(:,1:plotStep:end));
title('DC removed Raw Data')

upSampleData(Ch3DataObj);
nexttile
plot(Ch3DataObj.rawDataUpsampled(:,1:plotStep:end));
title('Upsampled Data')

alignWF_CFD(Ch3DataObj, 0.5)
removeSaturation(Ch3DataObj,dataLow,dataHigh);
avgData(Ch3DataObj, avg)
nexttile
plot(Ch3DataObj.averagedData(:,1:plotStep:end));
title('Averaged Data')

removeDCBG(Ch3DataObj,bgLow,bgHigh,0.03);
nexttile
plot(Ch3DataObj.preProcessedData(:,1:plotStep:end));
title('BG removed Data')

truncateData(Ch3DataObj,Truncation);
nexttile
plot(Ch3DataObj.dataT(:,1:plotStep:end));
title('Truncated Data')

runDeconLG(Ch3DataObj,exclude,LaguerreOrder,alpha);

runDeconExp(Ch3DataObj,expOrder,[],[],[],exclude);

%% save result
saveName = fullfile(filepath,[name '_12.5GS.mat']);
save(saveName,'dataInfoObj','Ch1DataObj','Ch2DataObj','Ch3DataObj')

close all

% %% plot fitting
% plotObj = Ch3DataObj;
% Idx = round(rand()*plotObj.numOfAvgWFs);
% rawDataToPlot = get(plotObj,'wf_aligned',Idx);
% rawDataMax = max(rawDataToPlot);
% fitToPlotLG = get(plotObj,'fit',Idx);
% fitToPlotExp = get(plotObj.expDeconObj,'fit',Idx);
% residueToPlotLG = rawDataToPlot-fitToPlotLG;
% residueToPlotExp = rawDataToPlot-fitToPlotExp;
% gain = plotObj.gain(Idx);
% CtrlV = plotObj.CtrlV(Idx);
% LT_LG =  plotObj.Lg_LTs(Idx);
% irfIdx = plotObj.irfIdx(Idx);
% irf = plotObj.APDObj.irfTNorm(:,irfIdx);
% t = plotObj.dtUp*(1:length(rawDataToPlot))';
% rawColor = [0 0 1 0.5];
% 
% figure
% plot(t,rawDataToPlot,'.-','Color', rawColor,'LineWidth', 1, 'MarkerSize',10)
% hold on
% plot(t(1:length(irf)),rawDataMax/max(irf)*irf,'Color',[0 0 0 0.5],'LineWidth', 1)
% plot(t,fitToPlotLG,'m','LineWidth', 1.2)
% plot(t,fitToPlotExp,'g-.','LineWidth', 1.2)
% hold off
% lgd = legend('Raw Data', 'iRF', 'LG fit','mExp fit');
% lgd.FontSize = 8;
% title_temp = sprintf('Point %2d, Laguerre Lifetime %3.2f, Gain %4.0f.',Idx,LT_LG,gain);
% % title_temp = [1 title_temp];
% title(title_temp)
% % ylim([-0.2 Ch1DataObj.dataHigh+0.1])
% 
% %%
% LifetimeOld = Ch3DataObj.Lg_LTs;
% Lifetime = Ch3DataObj_old.Lg_LTs;
% figure
% tiledlayout(1,2)
% nexttile
% scatter(Lifetime,Lifetime-LifetimeOld)
% nexttile
% histogram(Lifetime-LifetimeOld)
end