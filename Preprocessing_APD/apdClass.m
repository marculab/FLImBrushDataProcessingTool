classdef apdClass < handle
    % APD class to store APD detector data
    properties (Access = public)
        filePath = ''; % path to apd file
        serialNumber = ''; % apd serial number
        dateModified = ''; % last modified date
        user = ''; % user
        gainV % control voltage for gian, 1D vector
        apdGain % APD gain, 1D vector, must be the same size as gianV
        irfV % control voltage for irf, 1D vector of size M
        irf % irf matrix, 2D matrix of N by M
        irfTNorm % truncated and AUC normalized irf
        irfRawdt % time resolution of irf measurement
        irfdt % % time resolution of irf after interpolation
    end
    
    methods
        % contructor
        function apdObj = apdClass(pathToTdmsFile)
            apdObj.filePath = pathToTdmsFile;
        end
        
        function creatFromPath(apdObj)
            [~,~,ext] = fileparts(apdObj.filePath);
            switch ext
                case '.tdms'
                    loadTDMS(apdObj);
                case '.mat'
                    loadMat(apdObj);
            end
        end
        
        function loadMat(apdObj)
            temp = load(apdObj.filePath);
            apdObj.apdGain = temp.gain;
            apdObj.gainV = temp.gainV;
            apdObj.dateModified = temp.dateModified;
            apdObj.irf = temp.irf;
            apdObj.irfV = temp.irfV;
            apdObj.serialNumber = temp.serialNumber;
            apdObj.user = temp.user;
            apdObj.irfRawdt = temp.irfRawdt;
            apdObj.irfdt = temp.irfdt;
        end
        
        function loadTDMS(apdObj)
            [output,~] = TDMS_getStruct(apdObj.filePath);
            % set file name as serialNumber
            apdObj.serialNumber = output.Props.Serial_number;
            apdObj.dateModified = output.Props.DateTime;
            apdObj.user = output.Props.Author;
            % load in gain data
            gainVTemp =output.Gain.ControlV__V_.data;
            apdGainTemp = output.Gain.Gain.data;
            [gainVTemp,ia,~] = unique(gainVTemp);
            apdGainTemp = apdGainTemp(ia);
            apdObj.gainV = gainVTemp';
            apdObj.apdGain = apdGainTemp';
                       
            % load irf
            apdObj.irfRawdt = output.iRF.Props.dt;
            iRFName = fieldnames(output.iRFRaw);
            numOfIrf = length(iRFName)-2;
            try
            lengthOfIRF = output.iRFRaw.Props.WFLength;
            catch
               lengthOfIRF = output.iRF.Props.DataLength; 
            end
            try
            numOfiRFPerV = output.iRFRaw.Props.NumOfWFsPerV;
            catch
                numOfiRFPerV = output.iRFRaw.Props.DataLength/lengthOfIRF;
            end
            iRF = zeros(lengthOfIRF,numOfIrf);
            iRFV = zeros(numOfIrf,1);
            
            find_c = strfind(iRFName, 'c');
            iRFName(cellfun('isempty', find_c)) = [];
            iRFNameNumOnly = erase(iRFName,'_');
            iRFRawtemp = zeros(lengthOfIRF,numOfiRFPerV);
            for i = 1:numOfIrf
                iRFV(i) = sscanf(iRFNameNumOnly{i},'c%d')*0.01;
                iRFTemp = getfield(output.iRFRaw,iRFName{i});
                iRFRawtemp = reshape(iRFTemp.data,lengthOfIRF,numOfiRFPerV);
%                 figure; plot(iRFRawtemp);title('Before alignment');xlim([200 230])
                iRFRawtemp = alignWaveform_CFDNew(iRFRawtemp, 5, 0.2,0.5);
%                 figure; plot(iRFRawtemp);title('After alignment');xlim([200 230])
%                 for m = 1:numOfiRFPerV
%                     iRFRawtempUpSampled(:,m) = interp(iRFRawtemp(:,m),2);
%                 end
%                 figure;plot(iRFRawtempUpSampled);xlim([430 450])
%                 iRFRawtempUpSampled = alignWaveform_CFDNew(iRFRawtempUpSampled, 5, apdObj.irfdt, 0.5); % CFD alignment    
%                 figure;plot(iRFRawtempUpSampled);xlim([430 450])
%                 outlierIdx = isoutlier(max(iRFRawtemp));
                iRF(:,i) = mean(iRFRawtemp,2);
            end
            
            [iRFV,I] = sort(iRFV,'ascend');
            iRF = iRF(:,I);
%             iRF2G = interp1(Ch2CtrlV, Ch2Gain,Ch2iRFV,'makima');
            DC = mean(iRF(1:100,:));
            iRF = iRF-DC;
            apdObj.irfV = iRFV;
            apdObj.irf = iRF;
        end
    end
end
