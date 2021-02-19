classdef FLImDataClass < handle
    
    properties (Access = public)
        filePath; % path to data file
        numOfPoints; % number of raw data points
        WFLength; %waveform length
        ch1RawWF;
        ch2RawWF;
        ch3RawWF;
        V1;
        V2;
        V3;
        dataAvg; % number of waveforms to average
    end
    
    methods (Access = public)
        function FLImDataObj = FLImDataClass(pathToTdmsFile)
            FLImDataObj.filePath = pathToTdmsFile;
        end
        
        function loadFromfile(FLImDataObj)
            [output,~] = TDMS_getStruct(FLImDataObj.filePath);
            FLImDataObj.numOfPoints = output.Channel1.Props.Num_of_points;
            FLImDataObj.WFLength = output.Channel1.Props.Data_Length;
            FLImDataObj.ch1RawWF = reshape(output.Channel1.Waveform.data,FLImDataObj.WFLength,FLImDataObj.numOfPoints);
            FLImDataObj.V1 = output.Channel1.Voltage.data; 
            FLImDataObj.ch1RawWF = circshift(FLImDataObj.ch1RawWF,-0,2); % circular shift to account for delay
            FLImDataObj.V1 = circshift(FLImDataObj.V1,-1); % circular shift to account for delay
            FLImDataObj.ch2RawWF = reshape(output.Channel2.Waveform.data,FLImDataObj.WFLength,FLImDataObj.numOfPoints);
            FLImDataObj.V2 = output.Channel2.Voltage.data;
            FLImDataObj.ch2RawWF = circshift(FLImDataObj.ch2RawWF,-0,2); % circular shift to account for delay
            FLImDataObj.V2 = circshift(FLImDataObj.V2,-1); % circular shift to account for delay
            FLImDataObj.ch3RawWF = reshape(output.Channel3.Waveform.data,FLImDataObj.WFLength,FLImDataObj.numOfPoints);
            FLImDataObj.V3 = output.Channel3.Voltage.data;
            FLImDataObj.ch3RawWF = circshift(FLImDataObj.ch3RawWF,-0,2); % circular shift to account for delay
            FLImDataObj.V3 = circshift(FLImDataObj.V3,-1); % circular shift to account for delay
            checkDataIntegraty(FLImDataObj)
        end
        
        function checkDataIntegraty(FLImDataObj)
            ch1Integraty = mod(FLImDataObj.numOfPoints,length(FLImDataObj.V1));
            ch2Integraty = mod(FLImDataObj.numOfPoints,length(FLImDataObj.V2));
            ch3Integraty = mod(FLImDataObj.numOfPoints,length(FLImDataObj.V3));
            if (ch1Integraty+ch2Integraty+ch3Integraty)>0
                warndlg('Data dimention missmatch, possible corrupted data','Warning');
            else
%                 msgbox('FLIm data loaded successfully!')
                FLImDataObj.dataAvg = FLImDataObj.numOfPoints/length(FLImDataObj.V1);
            end
            
        end
    end
end
