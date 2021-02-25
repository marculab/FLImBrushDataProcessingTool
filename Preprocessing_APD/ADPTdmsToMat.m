% covert APD detector tdms file to .mat file
%%
clear all
close all
clc
%% select tdms file
[file,path] = uigetfile('.tdms','MultiSelect','on');
numOfFiles = length(file);

for n = 1: numOfFiles
    [output,~] = TDMS_getStruct(fullfile(path, file{n}));
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
    iRFUpSampled = zeros(lengthOfIRF*2,numOfIrf);
    iRFV = zeros(numOfIrf,1);
    
    find_c = strfind(iRFName, 'c');
    iRFName(cellfun('isempty', find_c)) = [];
    iRFNameNumOnly = erase(iRFName,'_');
    iRFRawtemp = zeros(lengthOfIRF,numOfiRFPerV);
    iRFRawtempUpSampled = zeros(lengthOfIRF*2,numOfiRFPerV);
    irfdt = 0.5*irfRawdt;
    for i = 1:numOfIrf
        iRFV(i) = sscanf(iRFNameNumOnly{i},'c%d')*0.01;
        iRFTemp = getfield(output.iRFRaw,iRFName{i});
        iRFRawtemp = reshape(iRFTemp.data,lengthOfIRF,numOfiRFPerV);
        for m = 1:numOfiRFPerV % upsample rawa data
            iRFRawtempUpSampled(:,m) = interp(iRFRawtemp(:,m),2);
        end
        %                 figure; plot(iRFRawtemp);title('Before alignment');xlim([200 230])
        iRFRawtemp = alignWaveform_CFDNew(iRFRawtemp, 3, irfRawdt,0.5);
        iRF(:,i) = mean(iRFRawtemp,2);
        
        iRFRawtempUpSampled = alignWaveform_CFDNew(iRFRawtempUpSampled, 3, irfdt, 0.5); % CFD alignment
        %                 figure;plot(iRFRawtempUpSampled);xlim([430 450])
        %                 outlierIdx = isoutlier(max(iRFRawtemp));
        iRFUpSampled(:,i) = mean(iRFRawtempUpSampled,2);
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
    [~,saveName,~] = fileparts(file{n});
    saveName = [saveName '.mat'];
    save(fullfile(path,saveName),'serialNumber','dateModified','user','gainV','gain','irfV','irf','irfUpSampled','irfRawdt','irfdt');
end
%%
idx = 100;
figure
plot(irf(:,idx));
hold on
plot(irfUpSampled(:,idx));
hold off

