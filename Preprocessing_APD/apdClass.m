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
        irfNorm % AUC normalized irf
        irfT = []; % trauncated irf
    end
    
    methods
        % contructor
        function apdObj = apdClass(pathToTdmsFile)
            apdObj.filePath = pathToTdmsFile;
        end
        
        function creatFromPath(apdObj)
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
            iRFName = fieldnames(output.iRF);
            numOfIrf = length(iRFName)-2;
            temp = getfield(output.iRF,iRFName{end});
            lengthOfIRF = length(temp.data);
            iRF = zeros(lengthOfIRF,numOfIrf);
            iRFV = zeros(numOfIrf,1);
            
            find_c = strfind(iRFName, 'c');
            iRFName(cellfun('isempty', find_c)) = [];
            iRFNameNumOnly = erase(iRFName,'_');
            for i = 1:numOfIrf
                iRFV(i) = sscanf(iRFNameNumOnly{i},'c%d')*0.01;
                iRFTemp = getfield(output.iRF,iRFName{i});
                iRF(:,i) = iRFTemp.data;
            end
            
            [iRFV,I] = sort(iRFV,'ascend');
            iRF = iRF(:,I);
%             iRF2G = interp1(Ch2CtrlV, Ch2Gain,Ch2iRFV,'makima');
            DC = mean(iRF(1:100,:));
            iRF = iRF-DC;
            apdObj.irfV = iRFV;
            apdObj.irf = iRF;
            apdObj.irfNorm = iRF./sum(iRF);
            
        end
    end
end
