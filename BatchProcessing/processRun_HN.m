% code to batch process HN data

function processRun_HN(root, DataPath, BGPath, APD1Path, APD2Path, APD3Path, APD4Path, alpha1, alpha2, alpha3, alpha4, ChWidth, LGorder)
addpath(genpath('..\Algorithms'))
addpath(genpath('..\BatchProcessing'))


dataInfoObj = dataInfo('','','','','','','','','','','','','','');
[filepath,name,ext] = fileparts(APD1Path);
dataInfoObj.apd1Name = [name ext];
dataInfoObj.apd1Folder = filepath;
[filepath,name,ext] = fileparts(APD2Path);
dataInfoObj.apd2Name = [name ext];
dataInfoObj.apd2Folder = filepath;
[filepath,name,ext] = fileparts(APD3Path);
dataInfoObj.apd3Name = [name ext];
dataInfoObj.apd3Folder = filepath;
[filepath,name,ext] = fileparts(APD4Path);
dataInfoObj.apd4Name = [name ext];
dataInfoObj.apd4Folder = filepath;
[filepath,name,ext] = fileparts(DataPath);
dataInfoObj.dataFileNames = [name ext];
dataInfoObj.dataFolderName = filepath;
[filepath,name,ext] = fileparts(BGPath);
dataInfoObj.bgFileName = [name ext];
dataInfoObj.bgFileFolder = filepath;

Ch1Color = [121,0,141]/255; % channel 1 plot color
Ch2Color = [0,169,255]/255; % channel 2 plot color
Ch3Color = [129,255,0]/255; % channel 3 plot color
Ch4Color = [255,0,0]/255; % channel 4 plot color

%% set variables 
upSampleFactor = 5;
%bgLow = 540;
%bgHigh = 900;
DigitizerNoise = 0.0078; % ENOB = 8;
numOfWFtoPlot = 100;
avg = 4;
% ChWidth = 680;
laguerreOrder = LGorder;
% laguerreOrder = 12;
% alpha1 = 0.916;
% alpha2 = 0.916;
% alpha3 = 0.916;
% alpha4 = 0.965;
phasorOrder = 4;

%% load APD file
% APD1Path = 'E:\MyGitRepo\FLImBrushDataProcessingTool\APDDetectorFile\UCD_GBSF\M00570374.mat';
apd1Obj = apdClass(APD1Path);
creatFromPath(apd1Obj);

% APD2Path = 'E:\MyGitRepo\FLImBrushDataProcessingTool\APDDetectorFile\UCD_GBSF\M00570375.mat';
apd2Obj = apdClass(APD2Path);
creatFromPath(apd2Obj);

% APD3Path = 'E:\MyGitRepo\FLImBrushDataProcessingTool\APDDetectorFile\UCD_GBSF\M00570376.mat';
apd3Obj = apdClass(APD3Path);
creatFromPath(apd3Obj);

% APD4Path = 'E:\MyGitRepo\FLImBrushDataProcessingTool\APDDetectorFile\UCD_GBSF\M00570378.mat';
apd4Obj = apdClass(APD4Path);
creatFromPath(apd4Obj);

% load in BG
% [foldername_BG,filename_BG,ext] = fileparts(BGPath);
% filename_BG = [filename_BG ext];
bgObj = backGround(char(fullfile(root, BGPath)));
loadBG(bgObj,upSampleFactor);
bgObj.bgGain1 = interp1(apd1Obj.gainV,apd1Obj.apdGain,bgObj.CtrlV1,'spline');
bgObj.bgGain2 = interp1(apd2Obj.gainV,apd2Obj.apdGain,bgObj.CtrlV2,'spline');
bgObj.bgGain3 = interp1(apd3Obj.gainV,apd3Obj.apdGain,bgObj.CtrlV3,'spline');
bgObj.bgGain4 = interp1(apd4Obj.gainV,apd4Obj.apdGain,bgObj.CtrlV4,'spline');

fBG = figure;
tiledlayout(4,1)
nexttile
plot(bgObj.bgCh1,'Marker','.','LineStyle','--','Color', Ch1Color,'MarkerSize',8);
title(sprintf('Channel 1 Background, CtrlV = %.3f, Gain = %.3f',bgObj.CtrlV1, bgObj.bgGain1));
ylim([min(bgObj.bgCh1)-0.1 max(bgObj.bgCh1)+0.1])
nexttile
plot(bgObj.bgCh2,'Marker','.','LineStyle','--','Color', Ch2Color,'MarkerSize',8);
title(sprintf('Channel 2 Background, CtrlV = %.3f, Gain = %.3f',bgObj.CtrlV2, bgObj.bgGain2));
ylim([min(bgObj.bgCh2)-0.1 max(bgObj.bgCh2)+0.1])
nexttile
plot(bgObj.bgCh3,'Marker','.','LineStyle','--','Color', Ch3Color,'MarkerSize',8);
title(sprintf('Channel 3 Background, CtrlV = %.3f, Gain = %.3f',bgObj.CtrlV3, bgObj.bgGain3));
ylim([min(bgObj.bgCh3)-0.1 max(bgObj.bgCh3)+0.1])
nexttile
plot(bgObj.bgCh4,'Marker','.','LineStyle','--','Color', Ch4Color,'MarkerSize',8);
title(sprintf('Channel 4 Background, CtrlV = %.3f, Gain = %.3f',bgObj.CtrlV4, bgObj.bgGain4));
ylim([min(bgObj.bgCh4)-0.1 max(bgObj.bgCh4)+0.1])

%% load in raw data
FLImDataObj = FLImDataClass(char(fullfile(root, DataPath)));
%FLImDataObj = FLImDataClass(fullfile(root, DataPath));
loadFromfile(FLImDataObj)
Ch1DataObj = ChannelDataAPD(FLImDataObj.ch1RawWF, FLImDataObj.V1, FLImDataObj.dataAvg, apd1Obj, 0.4, bgObj.bgCh1, FLImDataObj.laserRepRate, upSampleFactor);
Ch2DataObj = ChannelDataAPD(FLImDataObj.ch2RawWF, FLImDataObj.V2, FLImDataObj.dataAvg, apd2Obj, 0.4, bgObj.bgCh2, FLImDataObj.laserRepRate, upSampleFactor);
Ch3DataObj = ChannelDataAPD(FLImDataObj.ch3RawWF, FLImDataObj.V3, FLImDataObj.dataAvg, apd3Obj, 0.4, bgObj.bgCh3, FLImDataObj.laserRepRate, upSampleFactor);
Ch4DataObj = ChannelDataAPD(FLImDataObj.ch4RawWF, FLImDataObj.V4, FLImDataObj.dataAvg, apd4Obj, 0.4, bgObj.bgCh4, FLImDataObj.laserRepRate, upSampleFactor);

%% remove DC and plot raw waveforms
removeDCData(Ch1DataObj);
removeDCData(Ch2DataObj);
removeDCData(Ch3DataObj);
removeDCData(Ch4DataObj);

fRaw = figure;
tiledlayout(4,1)

titleText1 = ['Ch1 raw FLIm data (no averaging) ' num2str(numOfWFtoPlot) ' waveforms'];
titleText2 = ['Ch2 raw FLIm data (no averaging) ' num2str(numOfWFtoPlot) ' waveforms'];
titleText3 = ['Ch3 raw FLIm data (no averaging) ' num2str(numOfWFtoPlot) ' waveforms'];
titleText4 = ['Ch4 raw FLIm data (no averaging) ' num2str(numOfWFtoPlot) ' waveforms'];

plotStep = size(Ch1DataObj.preProcessedData,2)/numOfWFtoPlot;
plotStep = floor(plotStep);

temp = sortData(Ch1DataObj, 'ascend');
tempMin = min(temp(:));
if isnan(tempMin)
    tempMin = 0;
end
nexttile
plot(temp(:,1:plotStep:end))
title(titleText1)
ylim([tempMin-0.1 2]);
xlim([0 size(temp,1)]);

temp = sortData(Ch2DataObj, 'ascend');
if isnan(tempMin)
    tempMin = 0;
end
nexttile
plot(temp(:,1:plotStep:end))
title(titleText2)
ylim([tempMin-0.1 2]);
xlim([0 size(temp,1)]);

temp = sortData(Ch3DataObj, 'ascend');
if isnan(tempMin)
    tempMin = 0;
end
nexttile
plot(temp(:,1:plotStep:end))
title(titleText3)
ylim([tempMin-0.1 2]);
xlim([0 size(temp,1)]);

temp = sortData(Ch4DataObj, 'ascend');
if isnan(tempMin)
    tempMin = 0;
end
nexttile
plot(temp(:,1:plotStep:end))
title(titleText4)
ylim([tempMin-0.1 2]);
xlim([0 size(temp,1)]);

%%---------------------------upsample and align data-------------------------------
% upsample data as defined by factor 5
upSampleData(Ch1DataObj);
upSampleData(Ch2DataObj);
upSampleData(Ch3DataObj);
upSampleData(Ch4DataObj);

% align waveforms using a CFD with threshold 5
alignWF_CFD(Ch1DataObj, 0.5, 180*upSampleFactor:size(Ch2DataObj.rawData,1)*upSampleFactor)
alignWF_CFD(Ch2DataObj, 0.5, 180*upSampleFactor:size(Ch2DataObj.rawData,1)*upSampleFactor)
alignWF_CFD(Ch3DataObj, 0.5, 180*upSampleFactor:size(Ch2DataObj.rawData,1)*upSampleFactor)
alignWF_CFD(Ch4DataObj, 0.5, 180*upSampleFactor:size(Ch2DataObj.rawData,1)*upSampleFactor)

% plot upsampled data
fUpSampled = figure;
tiledlayout(4,1)
% plotStep = size(Ch1DataObj.preProcessedData,2)/numOfWFtoPlot;
% plotStep = floor(plotStep);
temp = sortData(Ch1DataObj, 'ascend');
nexttile
plot(temp(:,1:plotStep:end));
ylim([-0.1 2])
title('Ch1 DC&BG removed upsampled waveforms (100 Waveform)')

temp = sortData(Ch2DataObj, 'ascend');
nexttile
plot(temp(:,1:plotStep:end));
ylim([-0.1 2])
title('Ch2 DC&BG removed upsampled waveforms (100 Waveform)')

temp = sortData(Ch3DataObj, 'ascend');
nexttile
plot(temp(:,1:plotStep:end));
ylim([-0.1 2])
title('Ch3 DC&BG removed upsampled waveforms (100 Waveform)')

temp = sortData(Ch4DataObj, 'ascend');
nexttile
plot(temp(:,1:plotStep:end));
ylim([-0.1 2])
title('Ch4 DC&BG removed upsampled waveforms (100 Waveform)')
drawnow

%% ---------------remove saturation or low data points---------------------
% Remove saturation from data, remove waveforms with peak higher than 1.4V
% and lower than 0.03V
[totalNumOfWf, goodNumOfWf, badNumOfWf] = removeSaturation(Ch1DataObj,0.03,1.4);
msg1 = sprintf('Total number of waveforms: %d.\nWaveforms left: %d.\nWaveforms removed: %d.',...
    totalNumOfWf,goodNumOfWf,badNumOfWf);
[totalNumOfWf, goodNumOfWf, badNumOfWf] = removeSaturation(Ch2DataObj,0.03,1.4);
msg2 = sprintf('Total number of waveforms: %d.\nWaveforms left: %d.\nWaveforms removed: %d.',...
    totalNumOfWf,goodNumOfWf,badNumOfWf);
[totalNumOfWf, goodNumOfWf, badNumOfWf] = removeSaturation(Ch3DataObj,0.03,1.4);
msg3 = sprintf('Total number of waveforms: %d.\nWaveforms left: %d.\nWaveforms removed: %d.',...
    totalNumOfWf,goodNumOfWf,badNumOfWf);
[totalNumOfWf, goodNumOfWf, badNumOfWf] = removeSaturation(Ch4DataObj,0.03,1.4);
msg4 = sprintf('Total number of waveforms: %d.\nWaveforms left: %d.\nWaveforms removed: %d.',...
    totalNumOfWf,goodNumOfWf,badNumOfWf);
ff = msgbox({'Channel 1:';msg1;'';'Channel 2:';msg2;'';'Channel 3:';msg3;'Channel 4:';msg4},'Data info','help');
drawnow
close(ff)

%% ----------------------------averge data---------------------------------
% average data by a factor of 4
avgData(Ch1DataObj, avg)
avgData(Ch2DataObj, avg)
avgData(Ch3DataObj, avg)
avgData(Ch4DataObj, avg)

%% -------------- ---------remove fiber background-------------------------
% remove fiber background, first estimate if fiber background 
% In Channel 1 only, the plateau of the background is checked if background 
% needs estimation, if estimated background is higher than replaces
% background

% Define range for background based on maximum waveform location
[~, maxPosition] = max(Ch1DataObj.averagedData);
avgmaxPosition = round(mean(nonzeros(maxPosition)));
bgLow = avgmaxPosition - 420;
bgHigh = avgmaxPosition - 50 ;

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
% nexttile
% plot(estimatedBG)
% hold on
% h=plot(Ch1DataObj.bg,'r-','LineWidth',2);
% h2=plot(estimatedBGAvg,'g-','LineWidth',2);
% ylim([-0.2 1.4])
% title('Estimated BG from raw data')
% legend([h, h2],'Measured BG', 'Estimated BG')

Ch1DataObj.bg_Estimated = estimatedBGAvg; % overwrite existing BG with estimated BG
%Ch1DataObj.BgEstimated = 1; % set BG estimation flag

removeDCBG(Ch1DataObj,bgLow,bgHigh,1, 0.03);
removeDCBG(Ch2DataObj,bgLow,bgHigh,2);
removeDCBG(Ch3DataObj,bgLow,bgHigh,3);
removeDCBG(Ch4DataObj,bgLow,bgHigh,4);

ChWidthInPoints = ChWidth/Ch1DataObj.dtUp;

fBGremove = figure;
tiledlayout(4,2)
%----------------------plot channel 1-----------------------------------
nexttile
plot(Ch1DataObj.averagedData(:,1:plotStep:end))
xline(bgLow,'r--','BG','LineWidth',2);
xline(bgHigh,'r--','BG','LineWidth',2);
yline(DigitizerNoise,'m--','LineWidth',2);
ylim([tempMin-0.1 2]);
title('Ch1 waveforms (100 Waveform)');

temp = sortData(Ch1DataObj, 'ascend');
tempMin = min(temp(:));
if isnan(tempMin)
    tempMin=0;
end
[~,maxIdx] = max(temp);
maxIdx = median(maxIdx);
nexttile
plot(temp(:,1:plotStep:end))
xline(maxIdx-68+1,'c--','Truncation','LineWidth',2);
xline(maxIdx-68+1+ChWidthInPoints,'c--','Truncation','LineWidth',2);
xline(bgLow,'r--','LineWidth',2);
xline(bgHigh,'r--','LineWidth',2);
yline(DigitizerNoise,'m--','LineWidth',2);
ylim([tempMin-0.1 2]);
title('Ch1 DC&BG removed waveforms (100 Waveform)');
%----------------------plot channel 2-----------------------------------
nexttile
plot(Ch2DataObj.averagedData(:,1:plotStep:end))
xline(bgLow,'r--','BG','LineWidth',2);
xline(bgHigh,'r--','BG','LineWidth',2);
yline(DigitizerNoise,'m--','LineWidth',2);
ylim([tempMin-0.1 2]);
title('Ch2 waveforms (100 Waveform)');

temp = sortData(Ch2DataObj, 'ascend');
tempMin = min(temp(:));
if isnan(tempMin)
    tempMin=0;
end
[~,maxIdx] = max(temp);
maxIdx = median(maxIdx);
nexttile
plot(temp(:,1:plotStep:end))
xline(maxIdx-68+1,'c--','Truncation','LineWidth',2);
xline(maxIdx-68+1+ChWidthInPoints,'c--','Truncation','LineWidth',2);
xline(bgLow,'r--','LineWidth',2);
xline(bgHigh,'r--','LineWidth',2);
yline(DigitizerNoise,'m--','LineWidth',2);
ylim([tempMin-0.1 2]);
title('Ch2 DC&BG removed waveforms (100 Waveform)');
%----------------------plot channel 3-----------------------------------
nexttile
plot(Ch3DataObj.averagedData(:,1:plotStep:end))
xline(bgLow,'r--','BG','LineWidth',2);
xline(bgHigh,'r--','BG','LineWidth',2);
yline(DigitizerNoise,'m--','LineWidth',2);
ylim([tempMin-0.1 2]);
title('Ch3 waveforms (100 Waveform)');

temp = sortData(Ch3DataObj, 'ascend');
tempMin = min(temp(:));
if isnan(tempMin)
    tempMin=0;
end
[~,maxIdx] = max(temp);
maxIdx = median(maxIdx);
nexttile
plot(temp(:,1:plotStep:end))
xline(maxIdx-68+1,'c--','Truncation','LineWidth',2);
xline(maxIdx-68+1+ChWidthInPoints,'c--','Truncation','LineWidth',2);
xline(bgLow,'r--','LineWidth',2);
xline(bgHigh,'r--','LineWidth',2);
yline(DigitizerNoise,'m--','LineWidth',2);
ylim([tempMin-0.1 2]);
title('Ch3 DC&BG removed waveforms (100 Waveform)');
%----------------------plot channel 4-----------------------------------
nexttile
plot(Ch4DataObj.averagedData(:,1:plotStep:end))
xline(bgLow,'r--','BG','LineWidth',2);
xline(bgHigh,'r--','BG','LineWidth',2);
yline(DigitizerNoise,'m--','LineWidth',2);
ylim([tempMin-0.1 2]);
title('Ch4 waveforms (100 Waveform)');

temp = sortData(Ch4DataObj, 'ascend');
tempMin = min(temp(:));
if isnan(tempMin)
    tempMin=0;
end
[~,maxIdx] = max(temp);
maxIdx = median(maxIdx);
nexttile
plot(temp(:,1:plotStep:end))
xline(maxIdx-68+1,'c--','Truncation','LineWidth',2);
xline(maxIdx-68+1+ChWidthInPoints,'c--','Truncation','LineWidth',2);
xline(bgLow,'r--','LineWidth',2);
xline(bgHigh,'r--','LineWidth',2);
yline(DigitizerNoise,'m--','LineWidth',2);
ylim([tempMin-0.1 2]);
title('Ch4 DC&BG removed waveforms (100 Waveform)');

%% -------------------------truncate data----------------------------------
% Truncate data to 54.5 ns = 680 sample points
truncateData(Ch1DataObj,ChWidth);
truncateData(Ch2DataObj,ChWidth);
truncateData(Ch3DataObj,ChWidth);
truncateData(Ch4DataObj,ChWidth);

fTruncation = figure;
tiledlayout(4,1)
nexttile
temp = sortTruncatedData(Ch1DataObj, 'ascend');
plot(temp(:,1:plotStep:end));
yline(DigitizerNoise,'m--','LineWidth',2);
ylim([-0.1 2])
title('Ch1 DC&BG removed truncated waveforms (100 Waveform)')
nexttile
temp = sortTruncatedData(Ch2DataObj, 'ascend');
plot(temp(:,1:plotStep:end));
yline(DigitizerNoise,'m--','LineWidth',2);
ylim([-0.1 2])
title('Ch2 DC&BG removed truncated waveforms (100 Waveform)')
nexttile
temp = sortTruncatedData(Ch3DataObj, 'ascend');
plot(temp(:,1:plotStep:end));
yline(DigitizerNoise,'m--','LineWidth',2);
ylim([-0.1 2])
title('Ch3 DC&BG removed truncated waveforms (100 Waveform)')
nexttile
temp = sortTruncatedData(Ch4DataObj, 'ascend');
plot(temp(:,1:plotStep:end));
yline(DigitizerNoise,'m--','LineWidth',2);
ylim([-0.1 2])
title('Ch4 DC&BG removed truncated waveforms (100 Waveform)')

%%  ---------------------Laguerre Decon -------------------------------
% run laguerre deconvolution with set LG order and alpha values

exclude = 425:550;
if any(exclude < 1)
    exclude = [];
end
f = waitbar(0,'Running Laguerre channel 1...');
runDeconLG(Ch1DataObj,exclude,laguerreOrder,alpha1);

waitbar(0.25,f,'Running Laguerre channel 2...');
runDeconLG(Ch2DataObj,exclude,laguerreOrder,alpha2);

waitbar(0.5,f,'Running Laguerre channel 3...');
runDeconLG(Ch3DataObj,exclude,laguerreOrder,alpha3);

waitbar(0.75,f,'Running Laguerre channel 3...');
runDeconLG(Ch4DataObj,exclude,laguerreOrder,alpha4);

%  --------------------- Phasor analysis -------------------------------
% run Phasor analysis for order of harmonic 4
%for H = 1:4 % run phasor for harmonic of 1 to 4
for H = 1:phasorOrder
    runPhasor(Ch1DataObj, H);
    runPhasor(Ch2DataObj, H);
    runPhasor(Ch3DataObj, H);
    runPhasor(Ch4DataObj, H);
end

CH1WFGainCorr = get(Ch1DataObj,'GainCorrectedWF');
CH2WFGainCorr = get(Ch2DataObj,'GainCorrectedWF');
CH3WFGainCorr = get(Ch3DataObj,'GainCorrectedWF');
CH4WFGainCorr = get(Ch4DataObj,'GainCorrectedWF');
%-------------------------------------------compute extended phasor-----------------------------------------------------
WFExtended = [CH1WFGainCorr;CH2WFGainCorr;CH3WFGainCorr;CH4WFGainCorr];
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
Spectrum(:,4) = sum(CH4WFGainCorr)';
SP = zeros(numOfWF,1);
for i = 1: numOfWF
    SP(i) = (ComputePhasor(Spectrum(i,:),0,1,0)+(1+1i))/2;
end
SP_G = real(SP);
SP_S = imag(SP);

waitbar(1,f,'Finished');

close(f)
delete(f)

%% ---------------------------plot fitting---------------------------------
Idx = round(Ch1DataObj.numOfAvgWFs/4);
% Idx = 268;
fFitting1 = plotDeconFitting(Ch1DataObj,'Channel 1',Idx);
fFitting2 = plotDeconFitting(Ch2DataObj,'Channel 2',Idx);
fFitting3 = plotDeconFitting(Ch3DataObj,'Channel 3',Idx);
fFitting4 = plotDeconFitting(Ch4DataObj,'Channel 4',Idx);

%% covert to single
Ch1DataObj.toSingle();
Ch2DataObj.toSingle();
Ch3DataObj.toSingle();
Ch4DataObj.toSingle();

%% ----------------------------saving data---------------------------------
[filepath,fileName,~] = fileparts(DataPath); 
saveFolderName = 'ALL_DECONVOLVED_FILES';
mkdir(fullfile(root,filepath,saveFolderName)); 

%DateString = char(datetime("today"));
%saveFileName = fullfile(filepath,fileName,DateString,'.mat');
saveFileName = strcat(fileName, '.mat');

saveFileFullPath = fullfile(root,filepath,saveFolderName,saveFileName);
save(saveFileFullPath, 'dataInfoObj','Ch1DataObj','Ch2DataObj','Ch3DataObj','Ch4DataObj','EOP_H1G','EOP_H1S','SP_G','SP_S','-v7.3');
saveDeconLite(Ch1DataObj,Ch2DataObj,Ch3DataObj,Ch4DataObj,saveFileFullPath)

saveFolderNameFig = 'ALL_DECONVOLVED_FILES\figures';
mkdir(fullfile(root,filepath,saveFolderNameFig)); 
saveFigFullPath = fullfile(root,filepath,saveFolderNameFig,fileName);
saveas(fBG,strcat(saveFigFullPath,'_Background'),'fig');
saveas(fRaw,strcat(saveFigFullPath,'_RawWaveFrom'),'fig');
saveas(fUpSampled,strcat(saveFigFullPath,'_UpSampledWF'),'fig');
saveas(fBGremove,strcat(saveFigFullPath,'_BGRemoval'),'fig');
saveas(fTruncation,strcat(saveFigFullPath,'_TruncatedData'),'fig');
saveas(fFitting1,strcat(saveFigFullPath,'_FittingCh1'),'fig');
saveas(fFitting2,strcat(saveFigFullPath,'_FittingCh2'),'fig');
saveas(fFitting3,strcat(saveFigFullPath,'_FittingCh3'),'fig');
saveas(fFitting4,strcat(saveFigFullPath,'_FittingCh4'),'fig');

plotAndSaveFig(Ch1DataObj, bgLow, bgHigh, fileName, saveFigFullPath, 'Ch1')
plotAndSaveFig(Ch2DataObj, bgLow, bgHigh, fileName, saveFigFullPath, 'Ch2')
plotAndSaveFig(Ch3DataObj, bgLow, bgHigh, fileName, saveFigFullPath, 'Ch3')
plotAndSaveFig(Ch4DataObj, bgLow, bgHigh, fileName, saveFigFullPath, 'Ch4')
close all
end