% code to run multiexpoential fit to FLImBRUSH output

%%
clear 
close all
clc

%% add path
% addpath(genpath('./MultiExp Generic')) %add path to exponential code
addpath(genpath('..'))

%% load file
root = 'D:\LoaclData\5ALAFLImTest\Subject032_20210811';
[DeConfile,DeConpath] = uigetfile([root '\*.mat'],'Please select DeCon file','MultiSelect','on');
if iscell(DeConfile)
    numOfFile = length(DeConfile);
else
    numOfFile = 1;
end
for j = 1:numOfFile
    if iscell(DeConfile)
        P = fullfile(DeConpath,DeConfile{j});
    else
        P = fullfile(DeConpath,DeConfile);
    end
    load(P)
    % [filepath,name,ext] = fileparts(P);
    % fileName =
    %% Channel 1, variabel tau
    numOfWFCh1 = Ch1DataObj.numOfAvgWFs; % get number of waveforms
    orderCh1 = 3;
    deconIdxCh1 = Ch1DataObj.deconIdx;
    numOfDeconWFCh1 =  length(deconIdxCh1);
    ATemp = zeros(numOfDeconWFCh1,orderCh1);
    TauTemp = zeros(numOfDeconWFCh1,orderCh1);
    LTexpTemp = zeros(numOfDeconWFCh1,1);
    specCh1 = get(Ch1DataObj,'wf_aligned');
    irfNormCh1 = Ch1DataObj.APDObj.irfTNorm;
    irfIdxCh1 = Ch1DataObj.irfIdx;
    irfIdxValueCh1 = unique(irfIdxCh1,'sorted');
    fitCh1Temp = zeros(size(specCh1));
    decayCh1Temp = zeros(size(specCh1));
    
    for i = 1:length(irfIdxValueCh1)
        idxTemp = irfIdxValueCh1(i);
        laser = irfNormCh1(:,idxTemp);
        laser = laser./sum(laser);
        idxSpecTemp = find(irfIdxCh1==idxTemp);
        specTempCh1 = specCh1(:,idxSpecTemp);
        [A, T, avglife, intensity, fitCh1Temp(:,idxSpecTemp), raw, decayCh1Temp(:,idxSpecTemp)] = multiexp_fit(specTempCh1,Ch1DataObj.dtUp,laser,orderCh1,[]);
        ATemp(idxSpecTemp,:) = A';
        TauTemp(idxSpecTemp,:) = T';
        LTexpTemp(idxSpecTemp) = avglife;
    end
    WTemp = ATemp./sum(ATemp,2);
    resCh1Temp = specCh1-fitCh1Temp;
    % WTemp = WTemp';
    %
    ACh1 = zeros(numOfWFCh1,orderCh1);
    TauCh1 = zeros(numOfWFCh1,orderCh1);
    WCh1 = zeros(numOfWFCh1,orderCh1);
    LTexpCh1 = zeros(numOfWFCh1,1);
    fitCh1 = zeros(size(fitCh1Temp,1),numOfWFCh1);
    decayCh1 = zeros(size(decayCh1Temp,1),numOfWFCh1);
    resCh1 = zeros(size(decayCh1Temp,1),numOfWFCh1);
    
    
    ACh1(deconIdxCh1,:) = ATemp;
    TauCh1(deconIdxCh1,:) = TauTemp;
    WCh1(deconIdxCh1,:) = WTemp;
    LTexpCh1(deconIdxCh1) = LTexpTemp;
    fitCh1(:,deconIdxCh1) = fitCh1Temp;
    decayCh1(:,deconIdxCh1) = decayCh1Temp;
    resCh1(:,deconIdxCh1) = resCh1Temp;
    
    ExpResultCh1.A = ACh1;
    ExpResultCh1.Tau = TauCh1;
    ExpResultCh1.LTexp = LTexpCh1;
    ExpResultCh1.W = WCh1;
    ExpResultCh1.Order = orderCh1;
    ExpResultCh1.Fit = fitCh1;
    ExpResultCh1.Decay = decayCh1;
    ExpResultCh1.res = resCh1;
    
    tausCh1Avg = mean(TauCh1)
    
    %% Channel 1, fixed tau
    numOfWFCh1 = Ch1DataObj.numOfAvgWFs; % get number of waveforms
    deconIdxCh1 = Ch1DataObj.deconIdx;
    numOfDeconWFCh1 =  length(deconIdxCh1);
    specCh1 = get(Ch1DataObj,'wf_aligned');
    irfNormCh1 = Ch1DataObj.APDObj.irfTNorm;
    irfIdxCh1 = Ch1DataObj.irfIdx;
    irfIdxValueCh1 = unique(irfIdxCh1,'sorted');
    
    orderCh1 = 4;
    
    ATempF = zeros(numOfDeconWFCh1,orderCh1);
    TauTempF = zeros(numOfDeconWFCh1,orderCh1);
    LTexpTempF = zeros(numOfDeconWFCh1,1);
    fitCh1TempF = zeros(size(specCh1));
    decayCh1TempF = zeros(size(specCh1));
    
    for i = 1:length(irfIdxValueCh1)
        idxTemp = irfIdxValueCh1(i);
        laser = irfNormCh1(:,idxTemp);
        laser = laser./sum(laser);
        idxSpecTemp = find(irfIdxCh1==idxTemp);
        specTempCh1 = specCh1(:,idxSpecTemp);
        [A, T, avglife, intensity, fitCh1TempF(:,idxSpecTemp), raw, decayCh1TempF(:,idxSpecTemp)] = multiexp_fit(specTempCh1,Ch1DataObj.dtUp,laser,orderCh1,[]);
        ATempF(idxSpecTemp,:) = A';
        TauTempF(idxSpecTemp,:) = T';
        LTexpTempF(idxSpecTemp) = avglife;
    end
    WTempF = ATempF./sum(ATempF,2);
    resCh1TempF = specCh1-fitCh1TempF;
    
    % WTemp = WTemp';
    %
    ACh1F = zeros(numOfWFCh1,orderCh1);
    TauCh1F = zeros(numOfWFCh1,orderCh1);
    WCh1F = zeros(numOfWFCh1,orderCh1);
    LTexpCh1F = zeros(numOfWFCh1,1);
    fitCh1F = zeros(size(fitCh1TempF,1),numOfWFCh1);
    decayCh1F = zeros(size(decayCh1TempF,1),numOfWFCh1);
    resCh1 = zeros(size(decayCh1Temp,1),numOfWFCh1);
    
    ACh1F(deconIdxCh1,:) = ATempF;
    TauCh1F(deconIdxCh1,:) = TauTempF;
    WCh1F(deconIdxCh1,:) = WTempF;
    LTexpCh1F(deconIdxCh1) = LTexpTempF;
    fitCh1F(:,deconIdxCh1) = fitCh1TempF;
    decayCh1F(:,deconIdxCh1) = decayCh1TempF;
    resCh1F(:,deconIdxCh1) = resCh1TempF;
    
    ExpResultCh1F.A = ACh1F;
    ExpResultCh1F.Tau = TauCh1F;
    ExpResultCh1F.LTexp = LTexpCh1F;
    ExpResultCh1F.W = WCh1F;
    ExpResultCh1F.Order = orderCh1;
    ExpResultCh1F.Fit = fitCh1F;
    ExpResultCh1F.Decay = decayCh1F;
    ExpResultCh1F.res = resCh1F;
    %% -----------------------plot result-------------------------------------
    plotIdxCh1 = round(rand*numOfDeconWFCh1);
    figure('Position',[100 50 1200 900]);
    tiledlayout(3,2)
    nexttile
    plot(specCh1(:,plotIdxCh1),'r.-')
    hold on
    plot(laser,'g')
    plot(fitCh1Temp(:,plotIdxCh1),'b')
    title(['Channel 1 Point ' num2str(plotIdxCh1) '/' num2str(numOfDeconWFCh1)])
    legend('Raw Data','iRF','Fitting')
    hold off
    grid on
    ylim([-0.1 1])
%     axis tight
    nexttile
    plot(specCh1(:,plotIdxCh1),'r.-')
    hold on
    plot(laser,'g')
    plot(fitCh1TempF(:,plotIdxCh1),'b')
    title(['Channel 1 Point ' num2str(plotIdxCh1) '/' num2str(numOfDeconWFCh1)])
    legend('Raw Data','iRF','Fitting')
    hold off
    grid on
    ylim([-0.1 1])
%     axis tight
    nexttile
    plot(decayCh1Temp(:,plotIdxCh1)./max(decayCh1Temp(:,plotIdxCh1)))
    title(sprintf('Lifetime: %.3f',LTexpTemp(plotIdxCh1)))
    grid on
    axis tight
    nexttile
    plot(decayCh1TempF(:,plotIdxCh1)./max(decayCh1TempF(:,plotIdxCh1)))
    title(sprintf('Lifetime: %.3f',LTexpTempF(plotIdxCh1)))
    grid on
    axis tight
    nexttile
    plot(resCh1Temp(:,plotIdxCh1))
    grid on
    axis tight
    nexttile
    plot(resCh1TempF(:,plotIdxCh1))
    grid on
    axis tight
    drawnow
    
    %% plot lifetime diff and residual diff
    LTFiff = -LTexpTemp+LTexpTempF;
    resDiff = -sum(resCh1Temp.^2)+sum(resCh1TempF.^2);
    scatter(resDiff,LTFiff)
    %% Channel 2
    numOfWFCh2 = Ch2DataObj.numOfAvgWFs; % get number of waveforms
    orderCh2 = 3;
    deconIdxCh2 = Ch2DataObj.deconIdx;
    numOfDeconWFCh2 =  length(deconIdxCh2);
    ATemp = zeros(numOfDeconWFCh2,orderCh2);
    TauTemp = zeros(numOfDeconWFCh2,orderCh2);
    LTexpTemp = zeros(numOfDeconWFCh2,1);
    specCh2 = get(Ch2DataObj,'wf_aligned');
    irfNormCh2 = Ch2DataObj.APDObj.irfTNorm;
    irfIdxCh2 = Ch2DataObj.irfIdx;
    irfIdxValueCh2 = unique(irfIdxCh2,'sorted');
    fitCh2Temp = zeros(size(specCh2));
    decayCh2Temp = zeros(size(specCh2));
    
    for i = 1:length(irfIdxValueCh2)
        idxTemp = irfIdxValueCh2(i);
        laser = irfNormCh2(:,idxTemp);
        laser = laser./sum(laser);
        idxSpecTemp = find(irfIdxCh2==idxTemp);
        specTempCh2 = specCh2(:,idxSpecTemp);
        [A, T, avglife, intensity, fitCh2Temp(:,idxSpecTemp), raw, decayCh2Temp(:,idxSpecTemp)] = multiexp_fit(specTempCh2,Ch2DataObj.dtUp,laser,orderCh2,[0.3211    1.9296    8.1976]);
        ATemp(idxSpecTemp,:) = A';
        TauTemp(idxSpecTemp,:) = T';
        LTexpTemp(idxSpecTemp) = avglife;
    end
    WTemp = ATemp./sum(ATemp,2);
    % WTemp = WTemp';
    %
    ACh2 = zeros(numOfWFCh2,orderCh2);
    TauCh2 = zeros(numOfWFCh2,orderCh2);
    WCh2 = zeros(numOfWFCh2,orderCh2);
    LTexpCh2 = zeros(numOfWFCh2,1);
    fitCh2 = zeros(size(fitCh2Temp,1),numOfWFCh2);
    decayCh2 = zeros(size(decayCh2Temp,1),numOfWFCh2);
    
    ACh2(deconIdxCh2,:) = ATemp;
    TauCh2(deconIdxCh2,:) = TauTemp;
    WCh2(deconIdxCh2,:) = WTemp;
    LTexpCh2(deconIdxCh2) = LTexpTemp;
    fitCh2(:,deconIdxCh2) = fitCh2Temp;
    decayCh2(:,deconIdxCh2) = decayCh2Temp;
    
    ExpResultCh2.A = ACh2;
    ExpResultCh2.Tau = TauCh2;
    ExpResultCh2.LTexp = LTexpCh2;
    ExpResultCh2.W = WCh2;
    ExpResultCh2.Order = orderCh2;
    ExpResultCh2.Fit = fitCh2;
    ExpResultCh2.Decay = decayCh2;
    
    %-----------------------plot result-------------------------------------
    plotIdxCh2 = round(rand*numOfDeconWFCh2);
    figure('Position',[200 200 1200 600]);
    tiledlayout(2,1)
    nexttile
    plot(specCh2(:,plotIdxCh2),'r.-')
    hold on
    plot(laser,'g')
    plot(fitCh2Temp(:,plotIdxCh2),'b')
    hold off
    legend('Raw Data','iRF','Fitting')
    title(['Channel 2 Point' num2str(plotIdxCh2) '/' num2str(numOfDeconWFCh2)])
    nexttile
    plot(decayCh2Temp(:,plotIdxCh2))
    drawnow
    
    %% Channel 3
    numOfWFCh3 = Ch3DataObj.numOfAvgWFs; % get number of waveforms
    orderCh3 = 3;
    
    deconIdxCh3 = Ch3DataObj.deconIdx;
    numOfDeconWFCh3 =  length(deconIdxCh3);
    ATemp = zeros(numOfDeconWFCh3,orderCh3);
    TauTemp = zeros(numOfDeconWFCh3,orderCh3);
    LTexpTemp = zeros(numOfDeconWFCh3,1);
    specCh3 = get(Ch3DataObj,'wf_aligned');
    irfNormCh3 = Ch3DataObj.APDObj.irfTNorm;
    irfIdxCh3 = Ch3DataObj.irfIdx;
    irfIdxValueCh3 = unique(irfIdxCh3,'sorted');
    fitCh3Temp = zeros(size(specCh3));
    decayCh3Temp = zeros(size(specCh3));
    
    for i = 1:length(irfIdxValueCh3)
        idxTemp = irfIdxValueCh3(i);
        laser = irfNormCh3(:,idxTemp);
        laser = laser./sum(laser);
        idxSpecTemp = find(irfIdxCh3==idxTemp);
        specTempCh3 = specCh3(:,idxSpecTemp);
        [A, T, avglife, intensity, fitCh3Temp(:,idxSpecTemp), raw, decayCh3Temp(:,idxSpecTemp)] = multiexp_fit(specTempCh3,Ch3DataObj.dtUp,laser,orderCh3,[0.5533    2.4558   12.1994]);
        ATemp(idxSpecTemp,:) = A';
        TauTemp(idxSpecTemp,:) = T';
        LTexpTemp(idxSpecTemp) = avglife;
    end
    WTemp = ATemp./sum(ATemp,2);
    % WTemp = WTemp';
    %
    ACh3 = zeros(numOfWFCh3,3);
    TauCh3 = zeros(numOfWFCh3,3);
    WCh3 = zeros(numOfWFCh3,3);
    LTexpCh3 = zeros(numOfWFCh3,1);
    fitCh3 = zeros(size(fitCh3Temp,1),numOfWFCh3);
    decayCh3 = zeros(size(decayCh3Temp,1),numOfWFCh3);
    
    ACh3(deconIdxCh3,:) = ATemp;
    TauCh3(deconIdxCh3,:) = TauTemp;
    WCh3(deconIdxCh3,:) = WTemp;
    LTexpCh3(deconIdxCh3) = LTexpTemp;
    fitCh3(:,deconIdxCh3) = fitCh3Temp;
    decayCh3(:,deconIdxCh3) = decayCh3Temp;
    
    
    ExpResultCh3.A = ACh3;
    ExpResultCh3.Tau = TauCh3;
    ExpResultCh3.LTexp = LTexpCh3;
    ExpResultCh3.W = WCh3;
    ExpResultCh3.Order = orderCh3;
    ExpResultCh3.Fit = fitCh3;
    ExpResultCh3.Decay = decayCh3;
    
    %-----------------------plot result-------------------------------------
    plotIdxCh3 = round(rand*numOfDeconWFCh3);
    figure('Position',[200 200 1200 600]);
    tiledlayout(2,1)
    nexttile
    plot(specCh3(:,plotIdxCh3),'r.-')
    hold on
    plot(laser,'g')
    plot(fitCh3Temp(:,plotIdxCh3),'b')
    hold off
    legend('Raw Data','iRF','Fitting')
    title(['Channel 3 Point' num2str(plotIdxCh3) '/' num2str(numOfDeconWFCh3)])
    nexttile
    plot(decayCh3Temp(:,plotIdxCh3))
    drawnow
    
    %% save data
    % answer = questdlg('Save data?', ...
    % 	'Save data?', ...
    % 	'Yes','No','No');
    answer = 'Yes';
    switch answer
        case 'Yes'
            save(P, 'calibrationObj','Ch1DataObj','Ch2DataObj','Ch3DataObj','dataInfoObj','ExpResultCh1','ExpResultCh2','ExpResultCh3')
            disp('Data saved successfully!')
        case 'No'
            disp('No data is saved!')
    end
    %%
    close all
end