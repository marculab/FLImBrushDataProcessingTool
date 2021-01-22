classdef AMDModel < handle
    % specification for channeldata and iIRF:
    % 1. channeldata was background & DC subtracted, truncated and gain
    % corrected
    % 2. iIRF is unit-integral scaled, must be longer than channel width. 
    % * 10-90 percent before-after peak positioning is suggested for both
    % channeldata and iIRF
    
    properties
        LTs %lifetimes
        INTs % intensities
        channeldata % all relevent parameters and data
        shift % laser shift amount
        channelWidth % channel width
        iRF % processed iRF 
    end
    methods
        % constructor
        function obj = AMDModel(channeldata_in)
        obj.channeldata = channeldata_in;
        obj.channeldata.data = obj.channeldata.data(1:400,:);
        obj.channelWidth = size(channeldata_in.data,1);
        
        if length(obj.channeldata.iIRF)<=obj.channelWidth
           warning('iRF is too short, AMD calculation might have large error.')
           obj.iRF = obj.channeldata.iIRF;
        else
        [~, I] = max(obj.channeldata.iIRF);
        points_front = floor(obj.channelWidth*0.1);
        points_back = obj.channelWidth - points_front-1;
        obj.iRF = obj.channeldata.iIRF(I-points_front:I+points_back);
        end
        
        end
        
        function calculateLT(obj)
           dataLT = meanDelay(obj.channeldata.data,obj.channeldata.dt);
           irfLT = meanDelay(obj.iRF,obj.channeldata.dt);
           obj.LTs = dataLT-irfLT;
           obj.INTs = sum(obj.channeldata.data);
        end
    end
    
end
