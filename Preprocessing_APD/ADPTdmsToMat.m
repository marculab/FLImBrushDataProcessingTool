% covert APD detector tdms file to .mat file
%%
clear all
close all
clc
%% select tdms file
[file,path] = uigetfile('.tdms','MultiSelect','on');

if iscell(file)
    numOfFiles = length(file);
else
    numOfFiles = 1;
end

upSampleFactor = 4;

for n = 1: numOfFiles
    if numOfFiles==1
        [output,~] = TDMS_getStruct(fullfile(path, file));
        [~,saveName,~] = fileparts(file);
    else
        [output,~] = TDMS_getStruct(fullfile(path, file{n}));
        [~,saveName,~] = fileparts(file{n});
    end
    serialNumber = output.Props.Serial_number;
    dateModified = output.Props.DateTime;
    user = output.Props.Author;
    % load in gain data
    gainVTemp =output.Gain.ControlV__V_.data;
    apdGainTemp = output.Gain.Gain.data;
    [gainVTemp,ia,~] = unique(gainVTemp);
    apdGainTemp = apdGainTemp(ia);
    gainV = gainVTemp';
    gain = apdGainTemp';
    
    % load irf
    irfRawdt = output.iRF.Props.dt;
    iRFName = fieldnames(output.iRFRaw);
    numOfIrf = length(iRFName)-2;
    lengthOfIRF = output.iRFRaw.Props.WFLength;
    numOfiRFPerV = output.iRFRaw.Props.NumOfWFsPerV;
    iRF = zeros(lengthOfIRF,numOfIrf);
    iRFUpSampled = zeros(lengthOfIRF*upSampleFactor,numOfIrf);
    iRFV = zeros(numOfIrf,1);
    
    find_c = strfind(iRFName, 'c');
    iRFName(cellfun('isempty', find_c)) = [];
    iRFNameNumOnly = erase(iRFName,'_');
    iRFRawtemp = zeros(lengthOfIRF,numOfiRFPerV);
    iRFRawtempUpSampled = zeros(lengthOfIRF*upSampleFactor,numOfiRFPerV);
    irfdt = irfRawdt/upSampleFactor;
    for i = 1:numOfIrf
        iRFV(i) = sscanf(iRFNameNumOnly{i},'c%d')*0.01;
        iRFTemp = getfield(output.iRFRaw,iRFName{i});
        iRFRawtemp = reshape(iRFTemp.data,lengthOfIRF,numOfiRFPerV);
        for m = 1:numOfiRFPerV % upsample rawa data
            iRFRawtempUpSampled(:,m) = interp(iRFRawtemp(:,m),upSampleFactor);
        end
        %                 figure; plot(iRFRawtemp);title('Before alignment');xlim([200 230])
        iRFRawtemp = alignWaveform_CFDNew(iRFRawtemp, 2.8, irfRawdt,0.5);
        iRF(:,i) = mean(iRFRawtemp,2);
                
        iRFRawtempUpSampledAligned = alignWaveform_CFDNew(iRFRawtempUpSampled, 2.8, irfdt, 0.5); % CFD alignment
        %         xlim([400 460])
                if i==225
                   figure;tiledlayout(1,2);nexttile;plot(iRFRawtempUpSampled);xlim([400 460]);
                   nexttile;plot(iRFRawtempUpSampledAligned);xlim([400 460]);
                end
        %                 outlierIdx = isoutlier(max(iRFRawtemp));
        iRFUpSampled(:,i) = mean(iRFRawtempUpSampledAligned,2);
    end
    
    [iRFV,I] = sort(iRFV,'ascend');
    iRF = iRF(:,I);
    iRFUpSampled = iRFUpSampled(:,I);
    %             iRF2G = interp1(Ch2CtrlV, Ch2Gain,Ch2iRFV,'makima');
    DC = mean(iRF(1:100,:));
    iRF = iRF-DC;
    irfV = iRFV;
    irf = iRF;
    DCUp = mean(iRFUpSampled(1:200,:)); % compute DC
    iRFUpSampled = iRFUpSampled-DCUp; % remove DC
    irfUpSampled = iRFUpSampled;
    saveName = [saveName '.mat'];
    save(fullfile(path,saveName),'serialNumber','dateModified','user','gainV','gain','irfV','irf','irfUpSampled','irfRawdt','irfdt');
end
%% plot result
figure
idx = 100;
x = 0:size(irf,1)-1;
plot(x,irf(:,idx));
hold on
xx = 0:1/upSampleFactor:size(irf,1)-1/upSampleFactor;
plot(xx,irfUpSampled(:,idx));
hold off
xlim([100 600])

