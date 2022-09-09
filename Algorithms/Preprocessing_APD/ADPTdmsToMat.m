% covert APD detector tdms file to .mat file
%%
clear 
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
    try 
    lengthOfIRF = output.iRFRaw.Props.WFLength;
    catch
        lengthOfIRF = 800;
    end
    
    try
    numOfiRFPerV = output.iRFRaw.Props.NumOfWFsPerV;
    catch
        numOfiRFPerV = output.iRFRaw.Props.DataLength/lengthOfIRF;
    end
    iRF = zeros(lengthOfIRF,numOfIrf);
    iRFUpSampled = zeros(lengthOfIRF*upSampleFactor,numOfIrf);
    iRFV = zeros(numOfIrf,1);
    
    find_c = strfind(iRFName, 'c');
    iRFName(cellfun('isempty', find_c)) = [];
    iRFNameNumOnly = erase(iRFName,'_');
    iRFRawtemp = zeros(lengthOfIRF,numOfiRFPerV);
    
    irfdt = irfRawdt/upSampleFactor;
    for i = 1:numOfIrf
        iRFRawtempUpSampled = zeros(lengthOfIRF*upSampleFactor,numOfiRFPerV);
        iRFV(i) = sscanf(iRFNameNumOnly{i},'c%d')*0.01;
        iRFTemp = getfield(output.iRFRaw,iRFName{i});
%         iRFRawtemp = reshape(iRFTemp.data,lengthOfIRF,numOfiRFPerV);
        iRFRawtemp = reshape(iRFTemp.data,lengthOfIRF,[]);
        iRFRawtemp(:,sum(iRFRawtemp)==0) = [];
        for m = 1:min(numOfiRFPerV,size(iRFRawtemp,2)) % upsample rawa data
            iRFRawtempUpSampled(:,m) = interp(iRFRawtemp(:,m),upSampleFactor);
        end
        %                 figure; plot(iRFRawtemp);title('Before alignment');xlim([200 230])
        iRFRawtemp = alignWaveform_CFDNew(iRFRawtemp, 3, irfRawdt,0.5);
        iRF(:,i) = mean(iRFRawtemp,2);
        iRFRawtempUpSampled(:,sum(iRFRawtempUpSampled)==0) = [];
        iRFRawtempUpSampledAligned = alignWaveform_CFDNew(iRFRawtempUpSampled, 3, irfdt, 0.5); % CFD alignment
        %         xlim([400 460])
%                 if i==225
%                    figure;tiledlayout(1,2);nexttile;plot(iRFRawtempUpSampled);
%                    nexttile;plot(iRFRawtempUpSampledAligned);
%                 end
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
%     save(fullfile(path,saveName),'serialNumber','dateModified','user','gainV','gain','irfV','irf','irfUpSampled','irfRawdt','irfdt');
end
f = msgbox('Data Processing Completed');
%% plot result
figure
idx = 170;
x = 0:size(irf,1)-1;
plot(x,irf(:,idx));
hold on
xx = 0:1/upSampleFactor:size(irf,1)-1/upSampleFactor;
plot(xx,irfUpSampled(:,idx));
hold off
% xlim([100 600])

