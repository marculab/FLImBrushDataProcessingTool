classdef sysCaliDataClass < handle
    % System calibration files class to store system calibration data
    
    properties (Access = public)
        filePath = ''; % path to apd file
        caliDateTime = ''; % last modified date time
        ch1V % control voltage channel 1
        ch2V % control voltage channel 2
        ch3V % control voltage channel 3
        ch1WF % channel 1 calibration waveform
        ch2WF % channel 2 calibration waveform
        ch3WF % channel 3 calibration waveform
    end
    
    methods
        % contructor
        function obj = sysCaliDataClass(pathToFile)
            obj.filePath = pathToFile;
        end
        
        function load(obj)
            [output,~] = TDMS_getStruct(obj.filePath);
            fileInfo = dir(obj.filePath);
            obj.caliDateTime = fileInfo.date;
            obj.ch1V = output.channel_1.Props.voltage;
            obj.ch2V = output.channel_2.Props.voltage;
            obj.ch3V = output.channel_3.Props.voltage;
            obj.ch1WF = output.channel_1.waveform.data';
            obj.ch2WF = output.channel_2.waveform.data';
            obj.ch3WF = output.channel_3.waveform.data';
        end
        
    end
    
end