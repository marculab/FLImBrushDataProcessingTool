classdef FLImDataClass < handle
    % Create object to read FLImBRUSH .tdms files
    
    properties (Access = public)
        filePath; % absolute or relative path to data .tdms file
        numOfPoints; % number of raw waveforms before averaging
        WFLength; % waveform length in points
        ch1RawWF; % channel 1 raw waveform data, 2D matrix
        ch2RawWF; % channel 2 raw waveform data, 2D matrix
        ch3RawWF; % channel 3 raw waveform data, 2D matrix
        ch4RawWF; % channel 4 raw waveform data, 2D matrix
        V1; % channel 1 raw control V, 1D column vector
        V2; % channel 2 raw control V, 1D column vector
        V3; % channel 3 raw control V, 1D column vector
        V4; % channel 3 raw control V, 1D column vector
        dataAvg; % number of waveforms to average used in LabView control software
        laserRepRate; % 355nm laser rep rate
    end
    
    methods (Access = public)
        function FLImDataObj = FLImDataClass(pathToTdmsFile)
            % constructor, create empty object with path to tdms file
            FLImDataObj.filePath = pathToTdmsFile;
        end
        
        function loadFromfile(FLImDataObj)
            % function to load data into the object from the tdms file 
            [output,~] = TDMS_getStruct(FLImDataObj.filePath);
            FLImDataObj.numOfPoints = output.Channel1.Props.Num_of_points;
            try
            FLImDataObj.dataAvg = output.Props.WFAverage;
            catch
                FLImDataObj.dataAvg = output.Channel1.Props.Num_of_points/length(output.Channel1.Voltage.data);
            end
            
            try
            FLImDataObj.laserRepRate = double(output.Props.LaserRepRate);
            catch
                prompt = 'No laser rep rate information available from data file. Please manually input laser rep rate!';
                dlgtitle = 'Please input laser rep rate';
                dims = [1 70];
                definput = {'250'};
                answer = inputdlg(prompt,dlgtitle,dims,definput);
                FLImDataObj.laserRepRate = str2double(answer{1});
            end
%             FLImDataObj.laserRepRate; % output laser reprate to 
            
            timeDelay = [4.2800;4.9563;5.0427;5.0427]; % time delay in ms
            shift = round(timeDelay./1000.*FLImDataObj.laserRepRate);
            FLImDataObj.WFLength = output.Channel1.Props.Data_Length;
            
            FLImDataObj.ch1RawWF = reshape(output.Channel1.Waveform.data,FLImDataObj.WFLength,FLImDataObj.numOfPoints);
            % check data and V dimention, if the same, load data, if not duplicate V
            
            if size(FLImDataObj.ch1RawWF,2)==length(output.Channel1.Voltage.data)
                FLImDataObj.V1 = output.Channel1.Voltage.data';
            else
                temp = repmat(output.Channel1.Voltage.data,FLImDataObj.dataAvg,1);
                FLImDataObj.V1 = temp(:);
            end
            FLImDataObj.ch1RawWF = circshift(FLImDataObj.ch1RawWF,-shift(1),2); % circular shift to account for delay
            %             FLImDataObj.V1 = circshift(FLImDataObj.V1,-1); % circular shift to account for delay
            FLImDataObj.ch2RawWF = reshape(output.Channel2.Waveform.data,FLImDataObj.WFLength,FLImDataObj.numOfPoints);
            if size(FLImDataObj.ch2RawWF,2)==length(output.Channel2.Voltage.data)
                FLImDataObj.V2 = output.Channel2.Voltage.data';
            else
                temp = repmat(output.Channel2.Voltage.data,FLImDataObj.dataAvg,1);
                FLImDataObj.V2 = temp(:);
            end
            FLImDataObj.ch2RawWF = circshift(FLImDataObj.ch2RawWF,-shift(2),2); % circular shift to account for delay
            %             FLImDataObj.V2 = circshift(FLImDataObj.V2,-1); % circular shift to account for delay
            FLImDataObj.ch3RawWF = reshape(output.Channel3.Waveform.data,FLImDataObj.WFLength,FLImDataObj.numOfPoints);
            FLImDataObj.V3 = output.Channel3.Voltage.data';
            if size(FLImDataObj.ch3RawWF,2)==length(output.Channel3.Voltage.data)
                FLImDataObj.V3 = output.Channel3.Voltage.data';
            else
                temp = repmat(output.Channel3.Voltage.data,FLImDataObj.dataAvg,1);
                FLImDataObj.V3 = temp(:);
            end
            FLImDataObj.ch3RawWF = circshift(FLImDataObj.ch3RawWF,-shift(3),2); % circular shift to account for delay
            %             FLImDataObj.V3 = circshift(FLImDataObj.V3,-1); % circular shift to account for delay
            if isfield(output,'Channel4')
            FLImDataObj.ch4RawWF = reshape(output.Channel4.Waveform.data,FLImDataObj.WFLength,FLImDataObj.numOfPoints);
            if size(FLImDataObj.ch4RawWF,2)==length(output.Channel4.Voltage.data)
                FLImDataObj.V4 = output.Channel4.Voltage.data';
            else
                temp = repmat(output.Channel4.Voltage.data,FLImDataObj.dataAvg,1);
                FLImDataObj.V4 = temp(:);
            end
            FLImDataObj.ch4RawWF = circshift(FLImDataObj.ch4RawWF,-shift(4),2); % circular shift to account for delay
            end
            checkDataIntegraty(FLImDataObj)
        end
        
        function checkDataIntegraty(FLImDataObj)
            if ne(length(FLImDataObj.V1),size(FLImDataObj.ch1RawWF,2))
                warndlg('Channel 1 Data and control V dimention missmatch, possible corrupted data','Warning');
            end
            if ne(length(FLImDataObj.V2),size(FLImDataObj.ch2RawWF,2))
                warndlg('Channel 2 Data and control V dimention missmatch, possible corrupted data','Warning');
            end
            if ne(length(FLImDataObj.V3),size(FLImDataObj.ch3RawWF,2))
                warndlg('Channel 3 Data and control V dimention missmatch, possible corrupted data','Warning');
            end
            if ~isempty(FLImDataObj.V4) && ne(length(FLImDataObj.V4),size(FLImDataObj.ch4RawWF,2))
                warndlg('Channel 4 Data and control V dimention missmatch, possible corrupted data','Warning');
            end
        end
    end
end
