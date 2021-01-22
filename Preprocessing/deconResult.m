classdef deconResult < handle
    
    properties (Access = private)
        rawDataObj = []; % raw data obj
        laguerreObjectCellArray = []; % 4 by 1 cell array
    end
    
    properties (Access = public)
        ltMap = cell(4,1);
        intMap = cell(4,1);
        snrMap = cell(4,1);
        ltMapInterp = cell(4,1);
        intMapInterp = cell(4,1);
        snrMapInterp = cell(4,1);
        timeStampMatrix = [];
    end
    
    methods
        % constructor
        function obj = deconResult(rawDataObjIn, laguerreObjectCellArrayIn)
            obj.rawDataObj = rawDataObjIn;
            obj.laguerreObjectCellArray = laguerreObjectCellArrayIn;
            
            % loop through all laguerre object and reconstruct images
            nLines = length(rawDataObjIn.lineStartIdx);
            nPointsPerLine = rawDataObjIn.lineStartIdx(2)-rawDataObjIn.lineStartIdx(1);
            % get time stamp ready for interplation
            obj.timeStampMatrix = reshape(rawDataObjIn.timeStamp,nPointsPerLine,nLines)';
            timeStampDelta = diff(obj.timeStampMatrix,1,2);
            timeStampDelta = mean(timeStampDelta(:));
            tSampling = (1:nPointsPerLine)*timeStampDelta; % creat grid for data points
            
            for i=1:4
                laguerreObjTemp = laguerreObjectCellArrayIn{i};
                obj.ltMap{i} = reshape(laguerreObjTemp.LTs,nPointsPerLine,nLines)';
                obj.intMap{i} = reshape(laguerreObjTemp.INTs,nPointsPerLine,nLines)';
                obj.snrMap{i} = reshape(laguerreObjTemp.channeldata.SNR,nPointsPerLine,nLines)';
                % use interp1 to interplate data to equally sampled point
                obj.ltMapInterp{i} = interpMatrix(obj.timeStampMatrix,obj.ltMap{i},tSampling);
                obj.intMapInterp{i} = interpMatrix(obj.timeStampMatrix,obj.intMap{i},tSampling);
                obj.snrMapInterp{i} = interpMatrix(obj.timeStampMatrix,obj.snrMap{i},tSampling);
            end
            
            
            
        end
        
    end
    
    
end


