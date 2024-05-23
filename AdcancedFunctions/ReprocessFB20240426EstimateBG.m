% code to batch reprocess FLImBRUSH .mat data
% Xiangnan 09/09/2022
function ReprocessFB20240426EstimateBG(input_mat_name,apd1Obj,apd2Obj,apd3Obj)
% input_mat_name must be absolute path

[filepath,name,~] = fileparts(input_mat_name);
runName = split(name,'_');
runName = join(runName(1:end-1),"_");
name = runName{1};
pathSplit = split(filepath,'\');
pathSplit{1}='D:';
filepath = join(pathSplit,'\');
filepath = filepath{1};
if ~exist(filepath, 'dir') % check folder exist, if not make folder
    mkdir(filepath)
end
%% load processed .mat file
temp = load(input_mat_name,'dataInfoObj','Ch1DataObj','Ch2DataObj','Ch3DataObj');
dataInfoObj = temp.dataInfoObj;
Ch1DataObj_old = temp.Ch1DataObj;
Ch2DataObj_old = temp.Ch2DataObj;
Ch3DataObj_old = temp.Ch3DataObj;

%% reload BG to new sampling rate
upSampleFactor = 5;
% bgObj = backGround(fullfile(bgFolder,dataInfoObj.bgFileName));
% loadBG(bgObj,upSampleFactor);
%% retrive common data to all channels
laserRepRateIn = Ch1DataObj_old.laserRepRate;
LVAvgIn = Ch1DataObj_old.dataAveragingLV;
% bgLow = 535;
% bgHigh = 900;
avg = Ch1DataObj_old.dataAveraging;
dataLow = 0.03;
dataHigh = 1.4;
ChWidthInPoints = 680;
Truncation = ChWidthInPoints*0.08; % same as V4
exclude = 450:550;
LaguerreOrder = 12;
alpha = 0.916;
plotStep = ceil(size(Ch1DataObj_old.rawData,2)/200);
DigitizerNoise = 0.0078; % ENOB = 8;

% expOrder = Ch1DataObj_old.expDeconObj.order;
%% retrive nessesery Ch1 data for reprocessing
Ch1rawData = Ch1DataObj_old.rawData;
Ch1CtrlV = Ch1DataObj_old.rawCtrlV;
Ch1APDObjIn = apd1Obj; % use new APD object
Ch1dtIn = Ch1DataObj_old.dtRaw;
bgCh1 = Ch1DataObj_old.bg;
% update channel1 data obj to new decon setting
Ch1DataObj = ChannelDataAPD(Ch1rawData, Ch1CtrlV, LVAvgIn, Ch1APDObjIn, Ch1dtIn, bgCh1, laserRepRateIn, upSampleFactor);
% nexttile
% plot(Ch1rawData(:,1:500:end))
% title('Raw Data')
f1=figure('Position',[50 50 1600 900]);
tiledlayout(4,2,'TileSpacing','tight') ;

removeDCData(Ch1DataObj);
nexttile
tempWF = sortData(Ch1DataObj, 'ascend');
plot(tempWF(:,1:plotStep:end));
ylim([-0.2 1.4])
title('Ch1 DC removed Raw Data')

nexttile
plot(Ch1DataObj.bg);
ylim([-0.2 1.4])
title('Ch1 measured BG')

upSampleData(Ch1DataObj);
% nexttile
% plot(Ch1DataObj.rawDataUpsampled(:,1:500:end));
% title('Upsampled Data')

%-----------------------estimate BG from data------------------------------
tempMax = max(Ch1DataObj.rawDataUpsampled);
tempIdx = find(tempMax>0.9&tempMax<1.1);
if isempty(tempIdx)
    [~,I] = sort(tempMax, 'descend');
    estimatedBG = Ch1DataObj.rawDataUpsampled(:,I(1:100));
else
    temp = Ch1DataObj.rawDataUpsampled(:,tempIdx);
    tempCtrlV = Ch1DataObj.rawCtrlV(tempIdx);
    [~,I] = sort(tempCtrlV,'descend');
    estimatedBG = temp(:,I(1:100));
end
estimatedBGAvg = mean(estimatedBG,2);
nexttile
plot(estimatedBG)
hold on
h=plot(Ch1DataObj.bg,'r-','LineWidth',2);
h2=plot(estimatedBGAvg,'g-','LineWidth',2);
ylim([-0.2 1.4])
title('Estimated BG from raw data')
legend([h, h2],'Measured BG', 'Estimated BG')

Ch1DataObj.bg = estimatedBGAvg; % overwrite existing BG with estimated BG
Ch1DataObj.BgEstimated = 1; % set BG estimation flag

alignWF_CFD(Ch1DataObj, 0.99, 170*upSampleFactor:300*upSampleFactor) % align using distal peak
nexttile
tempWF = sortData(Ch1DataObj, 'ascend');
plot(tempWF(:,1:plotStep:end));
ylim([-0.2 1.4])
title('Upsampled and CFD aligned Data')

removeSaturation(Ch1DataObj,dataLow,dataHigh);
avgData(Ch1DataObj, avg)

% BG remove based on maximum waveform location
[~, maxPosition] = max(Ch1DataObj.averagedData);
avgmaxPosition = round(mean(nonzeros(maxPosition)));
bgLow = avgmaxPosition - 420;
bgHigh = avgmaxPosition - 50 ;

nexttile
tempWF = sortData(Ch1DataObj, 'ascend');
plot(tempWF(:,1:plotStep:end));
xline(bgLow,'r--','BG','LineWidth',2);
xline(bgHigh,'r--','BG','LineWidth',2);
yline(DigitizerNoise,'m--','LineWidth',2);
ylim([-0.2 1.4])
title('Averaged Data')


removeDCBG(Ch1DataObj,bgLow,bgHigh,1);
nexttile
tempWF = sortData(Ch1DataObj, 'ascend');
[~,maxIdx] = max(tempWF);
maxIdx = median(maxIdx);
plot(tempWF(:,1:plotStep:end));
xline(maxIdx-68+1,'c--','Truncation','LineWidth',2);
xline(maxIdx+ChWidthInPoints-68,'c--','Truncation','LineWidth',2);
xline(bgLow,'r--','LineWidth',2);
xline(bgHigh,'r--','LineWidth',2);
yline(DigitizerNoise,'m--','LineWidth',2);
% ylim([-0.2 1.4])
title('BG removed Data')

truncateData(Ch1DataObj,Truncation);
nexttile
plot(Ch1DataObj.dataT(:,1:plotStep:end));
yline(DigitizerNoise,'m--','LineWidth',2);
title('Truncated Data, nagtive WF removed')
xline(min(exclude),'g--','LineWidth',2);
xline(max(exclude),'g--','LineWidth',2);

runDeconLG(Ch1DataObj,exclude,LaguerreOrder,alpha);
Ch1DataObj.stat_test = test_stats(get(Ch1DataObj,'wf_aligned'),get(Ch1DataObj,'fit'),Ch1DataObj.dtUp, 20);

nexttile
idxTemp = find(max(Ch1DataObj.dataT)>0.5&max(Ch1DataObj.dataT)<0.7);
if ~isempty(idxTemp)
idxTemp = idxTemp(round(length(idxTemp)/2));
else
    idxTemp = round(size(Ch1DataObj.averagedData,2)/2);
end
% idxTemp = 1387;
rawDataToPlot = get(Ch1DataObj,'wf_aligned',idxTemp);
rawDataMax = max(rawDataToPlot);
fitToPlotLG = get(Ch1DataObj,'fit',idxTemp);
gain = Ch1DataObj.gain(idxTemp);
CtrlV = Ch1DataObj.CtrlV(idxTemp);
LT_LG =  Ch1DataObj.Lg_LTs(idxTemp);
irfIdx = Ch1DataObj.irfIdx(idxTemp);
irf = Ch1DataObj.APDObj.irfTNorm(:,irfIdx);

t = (1:length(rawDataToPlot));
rawColor = [0 0 0 0.5];
plot(t,rawDataToPlot,'.-','Color', rawColor,'LineWidth', 1, 'MarkerSize',10)
hold('on')
plot((1:length(irf)),rawDataMax/max(irf)*irf,'Color',[0 0 1 0.5],'LineWidth', 1)
plot(t,fitToPlotLG,'m','LineWidth', 1.2)
hold('off')
lgd = legend('Raw Data', 'iRF', 'LG fit');
lgd.FontSize = 8;
title_temp = sprintf('Point %2d, Laguerre Lifetime %3.2f, Gain %4.0f.',idxTemp,LT_LG,gain);
title_temp = ['Ch1' title_temp];
title(title_temp)
% ylim([-0.4 plotObj.dataHigh+0.1])

savefig(f1,fullfile(filepath,[name ' Ch1.fig']),'compact')
exportgraphics(f1,fullfile(filepath,[name ' Ch1.png']),'Resolution',600)
% runDeconExp(Ch1DataObj,expOrder,[],[],[],exclude);


%% retrive nessesery Ch2 data for reprocessing
Ch2rawData = Ch2DataObj_old.rawData;
Ch2CtrlV = Ch2DataObj_old.rawCtrlV;
Ch2APDObjIn = apd2Obj; % use new APD object
Ch2dtIn = Ch2DataObj_old.dtRaw;
bgCh2 = Ch2DataObj_old.bg;
% update channel1 data obj to new decon setting
Ch2DataObj = ChannelDataAPD(Ch2rawData, Ch2CtrlV, LVAvgIn, Ch2APDObjIn, Ch2dtIn, bgCh2, laserRepRateIn, upSampleFactor);
% nexttile
% plot(Ch2rawData(:,1:500:end))
% title('Raw Data')
f2=figure('Position',[50 50 1600 900]);
tiledlayout(4,2,'TileSpacing','tight') ;

removeDCData(Ch2DataObj);
nexttile
tempWF = sortData(Ch2DataObj, 'ascend');
plot(tempWF(:,1:plotStep:end));
ylim([-0.2 1.4])
title('Ch2 DC removed Raw Data')

nexttile
plot(Ch2DataObj.bg);
ylim([-0.2 1.4])
title('Ch2 measured BG')

upSampleData(Ch2DataObj);
% nexttile
% plot(Ch2DataObj.rawDataUpsampled(:,1:500:end));
% title('Upsampled Data')

%-----------------------estimate BG from data------------------------------
tempMax = max(Ch2DataObj.rawDataUpsampled);
tempIdx = find(tempMax>0.9&tempMax<1.1);
temp = Ch2DataObj.rawDataUpsampled(:,tempIdx);
tempCtrlV = Ch2DataObj.rawCtrlV(tempIdx);
[~,I] = sort(tempCtrlV,'descend');
estimatedBG = temp(:,I(1:100));

nexttile
plot(estimatedBG)
hold on
h=plot(Ch2DataObj.bg,'r-','LineWidth',2);
ylim([-0.2 1.4])
title('Estimated BG from raw data')
legend(h,'Measured BG')

alignWF_CFD(Ch2DataObj, 0.99, 170*upSampleFactor:300*upSampleFactor) % align using distal peak
nexttile
tempWF = sortData(Ch2DataObj, 'ascend');
plot(tempWF(:,1:plotStep:end));
ylim([-0.2 1.4])
title('Upsampled and CFD aligned Data')

removeSaturation(Ch2DataObj,dataLow,dataHigh);
avgData(Ch2DataObj, avg)
nexttile
tempWF = sortData(Ch2DataObj, 'ascend');
plot(tempWF(:,1:plotStep:end));
xline(bgLow,'r--','BG','LineWidth',2);
xline(bgHigh,'r--','BG','LineWidth',2);
yline(DigitizerNoise,'m--','LineWidth',2);
ylim([-0.2 1.4])
title('Averaged Data')

removeDCBG(Ch2DataObj,bgLow,bgHigh,1);
nexttile
tempWF = sortData(Ch2DataObj, 'ascend');
[~,maxIdx] = max(tempWF);
maxIdx = median(maxIdx);
plot(tempWF(:,1:plotStep:end));
xline(maxIdx-68+1,'c--','Truncation','LineWidth',2);
xline(maxIdx+ChWidthInPoints-68,'c--','Truncation','LineWidth',2);
xline(bgLow,'r--','LineWidth',2);
xline(bgHigh,'r--','LineWidth',2);
yline(DigitizerNoise,'m--','LineWidth',2);
% ylim([-0.2 1.4])
title('BG removed Data')

truncateData(Ch2DataObj,Truncation);
nexttile
plot(Ch2DataObj.dataT(:,1:plotStep:end));
yline(DigitizerNoise,'m--','LineWidth',2);
title('Truncated Data, nagtive WF removed')
xline(min(exclude),'g--','LineWidth',2);
xline(max(exclude),'g--','LineWidth',2);

runDeconLG(Ch2DataObj,exclude,LaguerreOrder,alpha);
Ch2DataObj.stat_test = test_stats(get(Ch2DataObj,'wf_aligned'),get(Ch2DataObj,'fit'),Ch2DataObj.dtUp, 20);

nexttile
% idxTemp = find(max(Ch2DataObj.dataT)>0.4&max(Ch2DataObj.dataT)<0.6);
% idxTemp = idxTemp(50);
% idxTemp = 1387;
rawDataToPlot = get(Ch2DataObj,'wf_aligned',idxTemp);
rawDataMax = max(rawDataToPlot);
fitToPlotLG = get(Ch2DataObj,'fit',idxTemp);
gain = Ch2DataObj.gain(idxTemp);
CtrlV = Ch2DataObj.CtrlV(idxTemp);
LT_LG =  Ch2DataObj.Lg_LTs(idxTemp);
irfIdx = Ch2DataObj.irfIdx(idxTemp);
irf = Ch2DataObj.APDObj.irfTNorm(:,irfIdx);

t = (1:length(rawDataToPlot));
rawColor = [0 0 0 0.5];
plot(t,rawDataToPlot,'.-','Color', rawColor,'LineWidth', 1, 'MarkerSize',10)
hold('on')
plot((1:length(irf)),rawDataMax/max(irf)*irf,'Color',[0 0 1 0.5],'LineWidth', 1)
plot(t,fitToPlotLG,'m','LineWidth', 1.2)
hold('off')
lgd = legend('Raw Data', 'iRF', 'LG fit');
lgd.FontSize = 8;
title_temp = sprintf('Point %2d, Laguerre Lifetime %3.2f, Gain %4.0f.',idxTemp,LT_LG,gain);
title_temp = ['Ch2' title_temp];
title(title_temp)
% ylim([-0.4 plotObj.dataHigh+0.1])

savefig(f2,fullfile(filepath,[name ' Ch2.fig']),'compact')
exportgraphics(f2,fullfile(filepath,[name ' Ch2.png']),'Resolution',600)

%% retrive nessesery Ch3 data for reprocessing
Ch3rawData = Ch3DataObj_old.rawData;
Ch3CtrlV = Ch3DataObj_old.rawCtrlV;
Ch3APDObjIn = apd3Obj; % use new APD object
Ch3dtIn = Ch3DataObj_old.dtRaw;
bgCh3 = Ch3DataObj_old.bg;
% update channel1 data obj to new decon setting
Ch3DataObj = ChannelDataAPD(Ch3rawData, Ch3CtrlV, LVAvgIn, Ch3APDObjIn, Ch3dtIn, bgCh3, laserRepRateIn, upSampleFactor);
% nexttile
% plot(Ch3rawData(:,1:500:end))
% title('Raw Data')
f3=figure('Position',[50 50 1600 900]);
tiledlayout(4,2,'TileSpacing','tight') ;

removeDCData(Ch3DataObj);
nexttile
tempWF = sortData(Ch3DataObj, 'ascend');
plot(tempWF(:,1:plotStep:end));
ylim([-0.2 1.4])
title('Ch3 DC removed Raw Data')

nexttile
plot(Ch3DataObj.bg);
ylim([-0.2 1.4])
title('Ch3 measured BG')

upSampleData(Ch3DataObj);
% nexttile
% plot(Ch3DataObj.rawDataUpsampled(:,1:500:end));
% title('Upsampled Data')

%-----------------------estimate BG from data------------------------------
tempMax = max(Ch3DataObj.rawDataUpsampled);
tempIdx = find(tempMax>0.9&tempMax<1.1);
temp = Ch3DataObj.rawDataUpsampled(:,tempIdx);
tempCtrlV = Ch3DataObj.rawCtrlV(tempIdx);
[~,I] = sort(tempCtrlV,'descend');
estimatedBG = temp(:,I(1:100));

nexttile
plot(estimatedBG)
hold on
h=plot(Ch3DataObj.bg,'r-','LineWidth',2);
ylim([-0.2 1.4])
title('Estimated BG from raw data')
legend(h,'Measured BG')

alignWF_CFD(Ch3DataObj, 0.99, 170*upSampleFactor:300*upSampleFactor) % align using distal peak
nexttile
tempWF = sortData(Ch3DataObj, 'ascend');
plot(tempWF(:,1:plotStep:end));
ylim([-0.2 1.4])
title('Upsampled and CFD aligned Data')

removeSaturation(Ch3DataObj,dataLow,dataHigh);
avgData(Ch3DataObj, avg)
nexttile
tempWF = sortData(Ch3DataObj, 'ascend');
plot(tempWF(:,1:plotStep:end));
xline(bgLow,'r--','BG','LineWidth',2);
xline(bgHigh,'r--','BG','LineWidth',2);
yline(DigitizerNoise,'m--','LineWidth',2);
ylim([-0.2 1.4])
title('Averaged Data')

removeDCBG(Ch3DataObj,bgLow,bgHigh,1);
nexttile
tempWF = sortData(Ch3DataObj, 'ascend');
[~,maxIdx] = max(tempWF);
maxIdx = median(maxIdx);
plot(tempWF(:,1:plotStep:end));
xline(maxIdx-68+1,'c--','Truncation','LineWidth',2);
xline(maxIdx+ChWidthInPoints-68,'c--','Truncation','LineWidth',2);
xline(bgLow,'r--','LineWidth',2);
xline(bgHigh,'r--','LineWidth',2);
yline(DigitizerNoise,'m--','LineWidth',2);
% ylim([-0.2 1.4])
title('BG removed Data')

truncateData(Ch3DataObj,Truncation);
nexttile
plot(Ch3DataObj.dataT(:,1:plotStep:end));
yline(DigitizerNoise,'m--','LineWidth',2);
title('Truncated Data, nagtive WF removed.')
xline(min(exclude),'g--','LineWidth',2);
xline(max(exclude),'g--','LineWidth',2);

runDeconLG(Ch3DataObj,exclude,LaguerreOrder,alpha);
Ch3DataObj.stat_test = test_stats(get(Ch3DataObj,'wf_aligned'),get(Ch3DataObj,'fit'),Ch3DataObj.dtUp, 20);

nexttile
% idxTemp = find(max(Ch3DataObj.dataT)>0.4&max(Ch3DataObj.dataT)<0.6);
% idxTemp = idxTemp(50);
% idxTemp = 1387;
rawDataToPlot = get(Ch3DataObj,'wf_aligned',idxTemp);
rawDataMax = max(rawDataToPlot);
fitToPlotLG = get(Ch3DataObj,'fit',idxTemp);
gain = Ch3DataObj.gain(idxTemp);
CtrlV = Ch3DataObj.CtrlV(idxTemp);
LT_LG =  Ch3DataObj.Lg_LTs(idxTemp);
irfIdx = Ch3DataObj.irfIdx(idxTemp);
irf = Ch3DataObj.APDObj.irfTNorm(:,irfIdx);

t = (1:length(rawDataToPlot));
rawColor = [0 0 0 0.5];
plot(t,rawDataToPlot,'.-','Color', rawColor,'LineWidth', 1, 'MarkerSize',10)
hold('on')
plot((1:length(irf)),rawDataMax/max(irf)*irf,'Color',[0 0 1 0.5],'LineWidth', 1)
plot(t,fitToPlotLG,'m','LineWidth', 1.2)
hold('off')
lgd = legend('Raw Data', 'iRF', 'LG fit');
lgd.FontSize = 8;
title_temp = sprintf('Point %2d, Laguerre Lifetime %3.2f, Gain %4.0f.',idxTemp,LT_LG,gain);
title_temp = ['Ch3' title_temp];
title(title_temp)
% ylim([-0.4 plotObj.dataHigh+0.1])

savefig(f3,fullfile(filepath,[name ' Ch3.fig']),'compact')
exportgraphics(f3,fullfile(filepath,[name ' Ch3.png']),'Resolution',600)

%% run phasor
for H = 1:4 % run phasor for harmonic of 1 to 4
    runPhasor(Ch1DataObj, H);
    runPhasor(Ch2DataObj, H);
    runPhasor(Ch3DataObj, H);
end

CH1WFGainCorr = get(Ch1DataObj,'GainCorrectedWF');
CH2WFGainCorr = get(Ch2DataObj,'GainCorrectedWF');
CH3WFGainCorr = get(Ch3DataObj,'GainCorrectedWF');
%-------------------------------------------compute extended phasor-----------------------------------------------------
WFExtended = [CH1WFGainCorr;CH2WFGainCorr;CH3WFGainCorr];
numOfWF = size(WFExtended,2);
EOP = zeros(numOfWF,1);
for i = 1:numOfWF
    EOP(i) = ComputePhasor(WFExtended(:,i),0,1,0);
end
EOP_H1G = real(EOP);
EOP_H1S = imag(EOP);
%-------------------------------------------compute spectral phasor--------------------------------------------
Spectrum = zeros(numOfWF,3);
Spectrum(:,1) = sum(CH1WFGainCorr)';
Spectrum(:,2) = sum(CH2WFGainCorr)';
Spectrum(:,3) = sum(CH3WFGainCorr)';
SP = zeros(numOfWF,1);
for i = 1: numOfWF
    SP(i) = (ComputePhasor(Spectrum(i,:),0,1,0)+(1+1i))/2;
end
SP_G = real(SP);
SP_S = imag(SP);

%% save result
Ch1DataObj.toSingle();
Ch2DataObj.toSingle();
Ch3DataObj.toSingle();
DateString = char(datetime("today"));
saveName = fullfile(filepath,[name '_12.5GS_' DateString '.mat']);
save(saveName,'dataInfoObj','Ch1DataObj','Ch2DataObj','Ch3DataObj','EOP_H1G','EOP_H1S','SP_G','SP_S','-v7.3')
saveDeconLite(Ch1DataObj,Ch2DataObj,Ch3DataObj,[],saveName)
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