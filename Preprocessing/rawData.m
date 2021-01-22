classdef rawData < handle
    %UNTITLED2 Summary of this class goes here
    %   class that handles preprocessing of raw data
    % functionality:
    % 1. Load in raw datya
    % 2. DC removal
    % 3. BG removal
    
    properties (Access = private)
        raw_data_protected = [];
        bg_protected = [];
    end
    
    properties (Access = public)
        fileNames = [];
        fileFolder = [];
        data = [];
        bg = [];
        irf = [];
        timeStamp = [];
        gain = [];
        pos = [];
        lineStartIdx = [];
        ID = [];
        DCIndex = [1 100];
        BGIndex = [1 100];
        BGThrehold = 0;
        peaks = [];
        noise = [];
        dataThrehold=[0.03 0.74];
        badDataIdx = [];
        goodDataIdx = [];
        scalingFactor = [];
    end
    
    methods
        function obj = rawData(fileNamesIn,fileFolderIn,bg_In,irf_In)
            %rawData Construct an instance of this class
            %   Detailed explanation goes here
            obj.fileNames = fileNamesIn;
            obj.fileFolder = fileFolderIn;
            obj.irf = irf_In;
            obj.bg_protected = bg_In;
            obj.bg = obj.bg_protected;
        end
        
        function obj = loadData(obj)
            %loadData to load in all binary data
            %   Detailed explanation goes here
            [obj.raw_data_protected,obj.gain,obj.timeStamp,~,obj.pos,obj.lineStartIdx,obj.ID] = LoadRawTRIPLEXData(obj.fileNames,{obj.fileFolder});
            % make a copy of protected data for public access and
            % processing
            obj.data = obj.raw_data_protected;
            obj.noise = std(obj.data(end-100:end,:))';
        end
        
        function obj = dataThresholding(obj,dataThreholdIn)
            obj.dataThrehold = dataThreholdIn;
            obj.badDataIdx = detectSaturationFunction(obj.data,dataThreholdIn(1),dataThreholdIn(2));
            fullIdx = 1:size(obj.raw_data_protected,2);
            obj.goodDataIdx =setdiff(fullIdx,obj.badDataIdx)';
        end
        
        function obj = DCRemove(obj,DcIndexIn)
            %loadData to load in all binary data
            %   Detailed explanation goes here
            obj.DCIndex = DcIndexIn;
            dcMatrix = obj.raw_data_protected(obj.DCIndex(1):obj.DCIndex(2),:); % dc matrix from all data point
            dc = mean(dcMatrix(:));% calculate average value
            obj.data = obj.raw_data_protected-dc;
            if ~isempty(obj.bg)
            obj.bg.bg = obj.bg_protected.bg-mean(obj.bg_protected.bg(obj.DCIndex(1):obj.DCIndex(2)));
            end
        end
        
        function obj = BGRemove(obj,BGIndexIn,BGThreholdIn)
            % remove fiber backgrounf from raw data
            %check whether there is available bg structure
            if ~isempty(obj.bg)
            % 1. Check whether there is observable background
            obj.BGIndex = BGIndexIn;
            obj.BGThrehold = BGThreholdIn;
            bgAvg = mean(obj.data(obj.BGIndex(1):obj.BGIndex(2)));
            if bgAvg>obj.BGThrehold
            dataAUC = sum(obj.data(obj.BGIndex(1):obj.BGIndex(2),:));
            bgAUC = sum(obj.bg.bg(obj.BGIndex(1):obj.BGIndex(2)));
            obj.scalingFactor = dataAUC./bgAUC;
            obj.scalingFactor = obj.scalingFactor'; % make it colume vector
            obj.data = obj.data-obj.bg.bg*obj.scalingFactor';
            else
                msgbox('Data BG too low. BG remove skipped!');
            end
            else
                msgbox('No BG data. BG remove skipped!');
            end
        end
        
        function obj = findPeaks(obj)
           % set noise to 0
           % use good data to find peaks, ignore bad data
           dataMean = mean(obj.data(:,obj.goodDataIdx),2);
           dataMean(dataMean<0.02*max(dataMean)) = 0;
           nonZeroIdx = find(dataMean); 
           [~,E] = discretize(nonZeroIdx,4);
           E(E<=0) = 1; % lowest index to be 1
           thisPeaks = zeros(4,1);
           for i = 1:4
               temp = dataMean;
               temp([1:E(i) E(i+1):end])=0;
               [~,idx] = max(temp);
               thisPeaks(i) = idx;
           end
           obj.peaks = thisPeaks;
        end
        
        function out = truncateData(obj, channelWidth) % truncate all data into 4 channels
            out = cell(4,1);
            for i = 1:4
                startIdx = obj.peaks(i)-channelWidth*0.1;
                startIdx = floor(startIdx);
                endIdx = startIdx+channelWidth-1;
                out{i} = obj.data(startIdx:endIdx,:);
            end
        end
    % end of method section    
    end
end

