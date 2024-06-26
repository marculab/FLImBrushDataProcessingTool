classdef sysCaliDataClass < handle
    % System calibration files class to store system calibration data
    
    properties (Access = public)
        filePath = ''; % path to apd file
        SystemID % system id to track system
        caliDateTime = ''; % last modified date time
        ch1V % control voltage channel 1
        ch2V % control voltage channel 2
        ch3V % control voltage channel 3
        ch1WF % channel 1 calibration waveform
        ch2WF % channel 2 calibration waveform
        ch3WF % channel 3 calibration waveform
        ch1Gain % channel 1 gain
        ch2Gain % channel 2 gain
        ch3Gain % channel 3 gain
    end
    
    methods
        function obj = sysCaliDataClass(pathToFile)
            % contructor, create object that contains system clibration data.
            obj.filePath = pathToFile;
        end
        
        function load(obj)
            % load function to load data from .tdms into object
            [output,~] = TDMS_getStruct(obj.filePath);
            fileInfo = dir(obj.filePath);
            obj.SystemID = output.Props.SystemID;
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