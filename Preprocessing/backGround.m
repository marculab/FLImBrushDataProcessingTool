classdef backGround < handle
    
    properties
        bgCh1; %Array n by 1
        bgCh2; %Array n by 1
        bgCh3; %Array n by 1
        pathToBg % full path to bg
    end
    
    methods
        function obj = backGround(pathToBgIn)
            obj.pathToBg = pathToBgIn;
        end
        
        function loadBG(obj)
            [output,~] = TDMS_getStruct(obj.pathToBg);
            bg1 = output.Channel_1.AvgBgData.data;
            bg2 = output.Channel_2.AvgBgData.data;
            bg3 = output.Channel_3.AvgBgData.data;
            obj.bgCh1 = bg1';
            obj.bgCh2 = bg2';
            obj.bgCh3 = bg3';                       
        end
    end
    
end
