% code to batch process HN data

function processRun(DataPath, BGPath, APD1Path, APD2Path, APD3Path, alpha1, alpha2, alpha3, ChWidth, saveRootPath)
addpath(genpath('..\Algorithms'))
addpath(genpath('..\BatchProcessing'))
dataInfoObj = dataInfo('','','','','','','','','','','','');
[filepath,name,ext] = fileparts(APD1Path);
dataInfoObj.apd1Name = [name ext];
dataInfoObj.apd1Folder = filepath;
[filepath,name,ext] = fileparts(APD2Path);
dataInfoObj.apd2Name = [name ext];
dataInfoObj.apd2Folder = filepath;
[filepath,name,ext] = fileparts(APD3Path);
dataInfoObj.apd3Name = [name ext];
dataInfoObj.apd3Folder = filepath;
[filepath,name,ext] = fileparts(DataPath);
dataInfoObj.dataFileNames = [name ext];
dataInfoObj.dataFolderName = filepath;
[filepath,name,ext] = fileparts(BGPath);
dataInfoObj.bgFileName = [name ext];
dataInfoObj.bgFileFolder = filepath;
% saveRootPath = 'C:\Users\Xiangnan\Desktop\BatchProcessTest';
% DataPath = 'Z:\V2 Sacramento Database\SP4-Intuitive Surgical\2023-03-21-Visit\2023_03_22_PEM\2023_03_22_PEM_50.tdms';
% BGPath = 'Z:\V2 Sacramento Database\SP4-Intuitive Surgical\2023-03-21-Visit\2023_03_21_Baseline_Tissue\Background\Background01.tdms';
Ch1Color = [121,0,141]/255; % channel 1 plot color
Ch2Color = [0,169,255]/255; % channel 2 plot color
Ch3Color = [129,255,0]/255; % channel 3 plot color
upSampleFactor = 5;
bgLow = 610;
bgHigh = 980;
DigitizerNoise = 0.0078; % ENOB = 8;
numOfWFtoPlot = 100;
laguerreOrder = 12;
% ChWidth = 154;
% alpha1 = 0.916;
% alpha2 = 0.916;
% alpha3 = 0.916;
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

% load in BG
[foldername_BG,filename_BG,ext] = fileparts(BGPath);
filename_BG = [filename_BG ext];
bgObj = backGround(fullfile(foldername_BG,filename_BG));
loadBG(bgObj,upSampleFactor);
bgObj.bgGain1 = interp1(apd1Obj.gainV,apd1Obj.apdGain,bgObj.CtrlV1,'spline');
bgObj.bgGain2 = interp1(apd2Obj.gainV,apd2Obj.apdGain,bgObj.CtrlV2,'spline');
bgObj.bgGain3 = interp1(apd3Obj.gainV,apd3Obj.apdGain,bgObj.CtrlV3,'spline');

fBG = figure;
tiledlayout(3,1)
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

%% load in raw data
FLImDataObj = FLImDataClass(DataPath);
loadFromfile(FLImDataObj)
Ch1DataObj = ChannelDataAPD(FLImDataObj.ch1RawWF, FLImDataObj.V1, FLImDataObj.dataAvg, apd1Obj, 0.4, bgObj.bgCh1, FLImDataObj.laserRepRate, upSampleFactor);
Ch2DataObj = ChannelDataAPD(FLImDataObj.ch2RawWF, FLImDataObj.V2, FLImDataObj.dataAvg, apd2Obj, 0.4, bgObj.bgCh2, FLImDataObj.laserRepRate, upSampleFactor);
Ch3DataObj = ChannelDataAPD(FLImDataObj.ch3RawWF, FLImDataObj.V3, FLImDataObj.dataAvg, apd3Obj, 0.4, bgObj.bgCh3, FLImDataObj.laserRepRate, upSampleFactor);
removeDCData(Ch1DataObj);
removeDCData(Ch2DataObj);
removeDCData(Ch3DataObj);


fRaw = figure;
tiledlayout(3,1)

titleText1 = ['Ch1 raw FLIm data (no averaging) ' num2str(numOfWFtoPlot) ' waveforms'];
titleText2 = ['Ch2 raw FLIm data (no averaging) ' num2str(numOfWFtoPlot) ' waveforms'];
titleText3 = ['Ch3 raw FLIm data (no averaging) ' num2str(numOfWFtoPlot) ' waveforms'];

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

%---------------------------upsample data-------------------------------
upSampleData(Ch1DataObj);
upSampleData(Ch2DataObj);
upSampleData(Ch3DataObj);
alignWF_CFD(Ch1DataObj, 1, (80:180)*upSampleFactor)
alignWF_CFD(Ch2DataObj, 0.5, (180:size(Ch2DataObj.rawData,1))*upSampleFactor)
alignWF_CFD(Ch3DataObj, 0.5, (180:size(Ch2DataObj.rawData,1))*upSampleFactor)
% plot upsampled data
fUpSampled = figure;
tiledlayout(3,1)
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
drawnow
%% ---------------remove saturation or low data points---------------------
[totalNumOfWf, goodNumOfWf, badNumOfWf] = removeSaturation(Ch1DataObj,0.03,1.4);
msg1 = sprintf('Total number of waveforms: %d.\nWaveforms left: %d.\nWaveforms removed: %d.',...
    totalNumOfWf,goodNumOfWf,badNumOfWf);
[totalNumOfWf, goodNumOfWf, badNumOfWf] = removeSaturation(Ch2DataObj,0.03,1.4);
msg2 = sprintf('Total number of waveforms: %d.\nWaveforms left: %d.\nWaveforms removed: %d.',...
    totalNumOfWf,goodNumOfWf,badNumOfWf);
[totalNumOfWf, goodNumOfWf, badNumOfWf] = removeSaturation(Ch3DataObj,0.03,1.4);
msg3 = sprintf('Total number of waveforms: %d.\nWaveforms left: %d.\nWaveforms removed: %d.',...
    totalNumOfWf,goodNumOfWf,badNumOfWf);
ff = msgbox({'Channel 1:';msg1;'';'Channel 2:';msg2;'';'Channel 3:';msg3},'Data info','help');
drawnow
close(ff)
%% ----------------------------averge data---------------------------------
avg = 4;
avgData(Ch1DataObj, avg)
avgData(Ch2DataObj, avg)
avgData(Ch3DataObj, avg)

%% ---- redefine bgHigh and bgLow in relation to peak of decay ------------
[maxValue, maxPosition] = max(Ch1DataObj.averagedData);

avgmaxPosition = round(mean(nonzeros(maxPosition)));

bgLow = avgmaxPosition - 420;
bgHigh = avgmaxPosition - 50 ;

disp(['bgLow: ', num2str(bgLow), ', bgHigh: ', num2str(bgHigh)]);

%% -------------- ---------remove fiber background-------------------------
removeDCBG(Ch1DataObj,bgLow,bgHigh,0.03);
removeDCBG(Ch2DataObj,bgLow,bgHigh);
removeDCBG(Ch3DataObj,bgLow,bgHigh);

ChWidthInPoints = ChWidth/Ch1DataObj.dtUp;

fBGremove = figure;
tiledlayout(3,2)
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
xline(maxIdx-200+1,'c--','Truncation','LineWidth',2);
xline(maxIdx+ChWidthInPoints-200,'c--','Truncation','LineWidth',2);
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
xline(maxIdx-200+1,'c--','Truncation','LineWidth',2);
xline(maxIdx+ChWidthInPoints-200,'c--','Truncation','LineWidth',2);
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
xline(maxIdx-200+1,'c--','Truncation','LineWidth',2);
xline(maxIdx+ChWidthInPoints-200,'c--','Truncation','LineWidth',2);
xline(bgLow,'r--','LineWidth',2);
xline(bgHigh,'r--','LineWidth',2);
yline(DigitizerNoise,'m--','LineWidth',2);
ylim([tempMin-0.1 2]);
title('Ch3 DC&BG removed waveforms (100 Waveform)');

%% -------------------------truncate data----------------------------------
truncateData(Ch1DataObj,ChWidth);
truncateData(Ch2DataObj,ChWidth);
truncateData(Ch3DataObj,ChWidth);

fTruncation = figure;
tiledlayout(3,1)
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

%%  ---------------------Laguerre Decon -------------------------------
exclude = 550:650;
if any(exclude < 1)
    exclude = [];
end
f = waitbar(0,'Running Laguerre channel 1...');
runDeconLG(Ch1DataObj,exclude,laguerreOrder,alpha1);

waitbar(0.33,f,'Running Laguerre channel 2...');
runDeconLG(Ch2DataObj,exclude,laguerreOrder,alpha2);

waitbar(0.66,f,'Running Laguerre channel 3...');
runDeconLG(Ch3DataObj,exclude,laguerreOrder,alpha3);
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

waitbar(1,f,'Finished');

close(f)
delete(f)

%% ---------------------------plot fitting---------------------------------
Idx = round(Ch1DataObj.numOfAvgWFs/2);
% Idx = 268;
fFitting1 = plotDeconFitting(Ch1DataObj,'Channel 1',Idx);
fFitting2 = plotDeconFitting(Ch2DataObj,'Channel 2',Idx);
fFitting3 = plotDeconFitting(Ch3DataObj,'Channel 3',Idx);

%% ----------------------------saving data---------------------------------
[~,runName,~] = fileparts(DataPath); % get run name
saveFolderName = [runName '_DeCon'];
mkdir(saveRootPath,saveFolderName); % create folder with run name
saveFileName = [runName '.mat'];
saveFileFullPath = fullfile(saveRootPath,saveFolderName,saveFileName);
save(saveFileFullPath, 'dataInfoObj','Ch1DataObj','Ch2DataObj','Ch3DataObj','EOP_H1G','EOP_H1S','SP_G','SP_S','-v7.3');

saveDeconLite(Ch1DataObj,Ch2DataObj,Ch3DataObj,saveFileFullPath)

saveas(fBG,fullfile(saveRootPath,saveFolderName,'Background'),'fig');
saveas(fRaw,fullfile(saveRootPath,saveFolderName,'RawWaveFrom'),'fig');
saveas(fUpSampled,fullfile(saveRootPath,saveFolderName,'UpSampledWF'),'fig');
saveas(fBGremove,fullfile(saveRootPath,saveFolderName,'BGRemoval'),'fig');
saveas(fTruncation,fullfile(saveRootPath,saveFolderName,'TruncatedData'),'fig');
saveas(fFitting1,fullfile(saveRootPath,saveFolderName,'FittingCh1'),'fig');
saveas(fFitting2,fullfile(saveRootPath,saveFolderName,'FittingCh2'),'fig');
saveas(fFitting3,fullfile(saveRootPath,saveFolderName,'FittingCh3'),'fig');

close all
end