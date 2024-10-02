% code to batch reprocess FLImBRUSH .mat data
% Xiangnan 09/09/2022
function ReprocessFB20240822EstimateBG(input_mat_name,SaveHardDriveLetter,apd1Obj,apd2Obj,apd3Obj,apd4Obj)
% input_mat_name must be absolute path

[filepath,name,~] = fileparts(input_mat_name);
% runName = split(name,'_');
% runName = join(runName(1:end-1),"_");
% name = runName{1};
name = char(name);
pathSplit = split(filepath,'\');
pathSplit{1}=SaveHardDriveLetter;
filepath = join(pathSplit,'\');
filepath = filepath{1};
if ~exist(filepath, 'dir') % check folder exist, if not make folder
    mkdir(filepath)
end
%% load processed .mat file
temp = load(input_mat_name,'dataInfoObj','Ch1DataObj','Ch2DataObj','Ch3DataObj','Ch4DataObj');
if isfield(temp,'Ch4DataObj')
    numOfChannel = 4;
else
    numOfChannel = 3;
end
dataInfoObj = temp.dataInfoObj;
Ch1DataObj_old = temp.Ch1DataObj;
Ch2DataObj_old = temp.Ch2DataObj;
Ch3DataObj_old = temp.Ch3DataObj;
if numOfChannel==4
Ch4DataObj_old = temp.Ch4DataObj;
end
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
estimateBG(Ch1DataObj);

nexttile
h=plot(Ch1DataObj.bg,'r-','LineWidth',2);
hold on
h2=plot(Ch1DataObj.bg_Estimated,'g-','LineWidth',2);
ylim([-0.2 1.4])
title('Estimated BG from raw data')
legend([h, h2],'Measured BG', 'Estimated BG')

alignWF_CFD(Ch1DataObj, 0.8, 170*upSampleFactor:300*upSampleFactor) % align using distal peak
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
estimateBG(Ch2DataObj);

nexttile
h=plot(Ch2DataObj.bg,'r-','LineWidth',2);
hold on
h2=plot(Ch2DataObj.bg_Estimated,'g-','LineWidth',2);
ylim([-0.2 1.4])
title('Estimated BG from raw data')
legend([h, h2],'Measured BG', 'Estimated BG')

alignWF_CFD(Ch2DataObj, 0.5, 180*upSampleFactor:size(Ch2DataObj.rawData,1)*upSampleFactor) % align using distal peak
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

removeDCBG(Ch2DataObj,bgLow,bgHigh,2);
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

estimateBG(Ch3DataObj);

nexttile
h=plot(Ch3DataObj.bg,'r-','LineWidth',2);
hold on
h2=plot(Ch3DataObj.bg_Estimated,'g-','LineWidth',2);
ylim([-0.2 1.4])
title('Estimated BG from raw data')
legend([h, h2],'Measured BG', 'Estimated BG')

alignWF_CFD(Ch3DataObj, 0.5, 180*upSampleFactor:size(Ch3DataObj.rawData,1)*upSampleFactor) % align using distal peak
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

removeDCBG(Ch3DataObj,bgLow,bgHigh,3);
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

%% if channel 4 exist
if numOfChannel==4
    % retrive nessesery Ch3 data for reprocessing
    Ch4rawData = Ch4DataObj_old.rawData;
    Ch4CtrlV = Ch4DataObj_old.rawCtrlV;
    Ch4APDObjIn = apd4Obj; % use new APD object
    Ch4dtIn = Ch4DataObj_old.dtRaw;
    bgCh4 = Ch4DataObj_old.bg;
    % update channel1 data obj to new decon setting
    Ch4DataObj = ChannelDataAPD(Ch4rawData, Ch4CtrlV, LVAvgIn, Ch4APDObjIn, Ch4dtIn, bgCh4, laserRepRateIn, upSampleFactor);
    % nexttile
    % plot(Ch3rawData(:,1:500:end))
    % title('Raw Data')
    f4=figure('Position',[50 50 1600 900]);
    tiledlayout(4,2,'TileSpacing','tight') ;

    removeDCData(Ch4DataObj);
    nexttile
    tempWF = sortData(Ch4DataObj, 'ascend');
    plot(tempWF(:,1:plotStep:end));
    ylim([-0.2 1.4])
    title('Ch3 DC removed Raw Data')

    nexttile
    plot(Ch4DataObj.bg);
    ylim([-0.2 1.4])
    title('Ch3 measured BG')

    upSampleData(Ch4DataObj);
    % nexttile
    % plot(Ch3DataObj.rawDataUpsampled(:,1:500:end));
    % title('Upsampled Data')

    estimateBG(Ch4DataObj);

    nexttile
    h=plot(Ch4DataObj.bg,'r-','LineWidth',2);
    hold on
    h2=plot(Ch4DataObj.bg_Estimated,'g-','LineWidth',2);
    ylim([-0.2 1.4])
    title('Estimated BG from raw data')
    legend([h, h2],'Measured BG', 'Estimated BG')

    alignWF_CFD(Ch4DataObj, 0.5, 180*upSampleFactor:size(Ch4DataObj.rawData,1)*upSampleFactor) % align using distal peak
    nexttile
    tempWF = sortData(Ch4DataObj, 'ascend');
    plot(tempWF(:,1:plotStep:end));
    ylim([-0.2 1.4])
    title('Upsampled and CFD aligned Data')

    removeSaturation(Ch4DataObj,dataLow,dataHigh);
    avgData(Ch4DataObj, avg)
    nexttile
    tempWF = sortData(Ch4DataObj, 'ascend');
    plot(tempWF(:,1:plotStep:end));
    xline(bgLow,'r--','BG','LineWidth',2);
    xline(bgHigh,'r--','BG','LineWidth',2);
    yline(DigitizerNoise,'m--','LineWidth',2);
    ylim([-0.2 1.4])
    title('Averaged Data')

    removeDCBG(Ch4DataObj,bgLow,bgHigh,3);
    nexttile
    tempWF = sortData(Ch4DataObj, 'ascend');
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

    truncateData(Ch4DataObj,Truncation);
    nexttile
    plot(Ch4DataObj.dataT(:,1:plotStep:end));
    yline(DigitizerNoise,'m--','LineWidth',2);
    title('Truncated Data, nagtive WF removed.')
    xline(min(exclude),'g--','LineWidth',2);
    xline(max(exclude),'g--','LineWidth',2);

    runDeconLG(Ch4DataObj,exclude,LaguerreOrder,alpha);

    nexttile
    % idxTemp = find(max(Ch3DataObj.dataT)>0.4&max(Ch3DataObj.dataT)<0.6);
    % idxTemp = idxTemp(50);
    % idxTemp = 1387;
    rawDataToPlot = get(Ch4DataObj,'wf_aligned',idxTemp);
    rawDataMax = max(rawDataToPlot);
    fitToPlotLG = get(Ch4DataObj,'fit',idxTemp);
    gain = Ch4DataObj.gain(idxTemp);
    CtrlV = Ch4DataObj.CtrlV(idxTemp);
    LT_LG =  Ch4DataObj.Lg_LTs(idxTemp);
    irfIdx = Ch4DataObj.irfIdx(idxTemp);
    irf = Ch4DataObj.APDObj.irfTNorm(:,irfIdx);

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
    title_temp = ['Ch4' title_temp];
    title(title_temp)
    % ylim([-0.4 plotObj.dataHigh+0.1])

    savefig(f4,fullfile(filepath,[name ' Ch4.fig']),'compact')
    exportgraphics(f4,fullfile(filepath,[name ' Ch4.png']),'Resolution',600)
end
%% run phasor
for H = 1:4 % run phasor for harmonic of 1 to 4
    runPhasor(Ch1DataObj, H);
    runPhasor(Ch2DataObj, H);
    runPhasor(Ch3DataObj, H);
    if numOfChannel ==4
        runPhasor(Ch4DataObj, H);
    end
end

CH1WFGainCorr = get(Ch1DataObj,'GainCorrectedWF');
CH2WFGainCorr = get(Ch2DataObj,'GainCorrectedWF');
CH3WFGainCorr = get(Ch3DataObj,'GainCorrectedWF');
if numOfChannel ==4
    CH4WFGainCorr = get(Ch4DataObj,'GainCorrectedWF');
else
    CH4WFGainCorr = [];
end
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

%% save result
Ch1DataObj.toSingle();
Ch2DataObj.toSingle();
Ch3DataObj.toSingle();

DateString = char(datetime("31-Aug-2024"));
saveName = fullfile(filepath,[name '_' DateString '_12.5GS.mat']);
if numOfChannel==4
    Ch4DataObj.toSingle();
    save(saveName,'dataInfoObj','Ch1DataObj','Ch2DataObj','Ch3DataObj','Ch4DataObj','EOP_H1G','EOP_H1S','SP_G','SP_S','-v7.3')
    saveDeconLite(Ch1DataObj,Ch2DataObj,Ch3DataObj,Ch4DataObj,saveName)
else
    save(saveName,'dataInfoObj','Ch1DataObj','Ch2DataObj','Ch3DataObj','EOP_H1G','EOP_H1S','SP_G','SP_S','-v7.3')
    saveDeconLite(Ch1DataObj,Ch2DataObj,Ch3DataObj,[],saveName)
end

close all
end