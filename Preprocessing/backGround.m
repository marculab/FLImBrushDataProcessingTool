classdef backGround < handle
    
    properties
        bgCh1; %Array n by 1
        bgCh2; %Array n by 1
        bgCh3; %Array n by 1
        CtrlV1; % scaler
        CtrlV2; % scaler
        CtrlV3; % scaler
        bgGain1 = NaN; % scaler gain value
        bgGain2 = NaN; % scaler gain value
        bgGain3 = NaN; % scaler gain value
        pathToBg % full path to bg
    end
    
    methods
        function obj = backGround(pathToBgIn)
            obj.pathToBg = pathToBgIn;
        end
        
        function loadBG(obj)
            [output,~] = TDMS_getStruct(obj.pathToBg);
            try
            bg1 = output.Channel_1.AvgBgData.data;
            bg2 = output.Channel_2.AvgBgData.data;
            bg3 = output.Channel_3.AvgBgData.data;
            obj.bgCh1 = bg1';
            obj.bgCh2 = bg2';
            obj.bgCh3 = bg3';
            obj.CtrlV1 = output.Channel_1.CtrlV.data(1);
            obj.CtrlV2 = output.Channel_2.CtrlV.data(1);
            obj.CtrlV3 = output.Channel_3.CtrlV.data(1);
            catch
               bg1 = output.Background.Ch1.data;
               bg2 = output.Background.Ch2.data;
               bg3 = output.Background.Ch3.data;
               obj.bgCh1 = bg1';
               obj.bgCh2 = bg2';
               obj.bgCh3 = bg3';    
               obj.CtrlV1 = NaN;
               obj.CtrlV2 = NaN;
               obj.CtrlV3 = NaN;
            end
        end
    end
    
end
