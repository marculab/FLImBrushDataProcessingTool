classdef backGround < handle
    
    properties
        bgCh1; % BG upsampled
        bgCh1Raw; %Array n by 1
        bgCh1Aligned; %Array n by 1
        bgCh2; % BG upsampled
        bgCh2Raw; %Array n by 1
        bgCh2Aligned; %Array n by 1
        bgCh3; % BG upsampled
        bgCh3Raw; %Array n by 1
        bgCh3Aligned; %Array n by 1
        CtrlV1; % scaler
        CtrlV2; % scaler
        CtrlV3; % scaler
        bgGain1 = NaN; % scaler gain value
        bgGain2 = NaN; % scaler gain value
        bgGain3 = NaN; % scaler gain value
        pathToBg % full path to bg
        numOfWF; % number of waveforms
        wfLength; % length of waveform
        upSampleFactor = 2; % upsample factor
        dt = 0.4; % times resolution 0.4 ns/2.5 GS/s
    end
    
    methods
        function obj = backGround(pathToBgIn)
            obj.pathToBg = pathToBgIn;
        end
        
        function loadBG(obj)
            [output,~] = TDMS_getStruct(obj.pathToBg);
            try
                %---------------------------------------get raw waveforms-----------------------------------------------
                obj.numOfWF = length(output.Channel_1.CtrlV.data); % get number of waveforms
                obj.wfLength = length(output.Channel_1.AvgBgData.data); % get waveform length
                obj.bgCh1Raw = reshape(output.Channel_1.RawBgData.data,obj.wfLength,obj.numOfWF); % get raw waveforms
                obj.bgCh2Raw = reshape(output.Channel_2.RawBgData.data,obj.wfLength,obj.numOfWF); % get raw waveforms
                obj.bgCh3Raw = reshape(output.Channel_3.RawBgData.data,obj.wfLength,obj.numOfWF); % get raw waveforms
                %---------------------------------------upsample data-------------------------------------------------
                bgCh1Temp = zeros(obj.wfLength*obj.upSampleFactor,obj.numOfWF);
                bgCh2Temp = zeros(obj.wfLength*obj.upSampleFactor,obj.numOfWF);
                bgCh3Temp = zeros(obj.wfLength*obj.upSampleFactor,obj.numOfWF);
                for i = 1:obj.numOfWF
                    bgCh1Temp(:,i) = interp(obj.bgCh1Raw(:,i),obj.upSampleFactor);
                    bgCh2Temp(:,i) = interp(obj.bgCh2Raw(:,i),obj.upSampleFactor);
                    bgCh3Temp(:,i) = interp(obj.bgCh3Raw(:,i),obj.upSampleFactor);
                end
                obj.bgCh1Aligned = alignWaveform_CFDNew(bgCh1Temp-0.5, 2.8, obj.dt/obj.upSampleFactor,1);
                obj.bgCh2Aligned = alignWaveform_CFDNew(bgCh2Temp, 2.6, obj.dt/obj.upSampleFactor,0.6);
                obj.bgCh3Aligned = alignWaveform_CFDNew(bgCh3Temp, 2.8, obj.dt/obj.upSampleFactor,0.6);
                %--------------------------------------plot alignment result-------------------------------------------
%                 figure;
%                 tiledlayout(1,2)
%                 nexttile
%                 plot(bgCh1Temp)
%                 xlim([450 480])
%                 nexttile
%                 plot(obj.bgCh1Aligned)
%                 xlim([450 480])
%                 
%                 figure;
%                 tiledlayout(1,2)
%                 nexttile
%                 plot(bgCh2Temp)
%                 xlim([450 480])
%                 nexttile
%                 plot(obj.bgCh2Aligned)
%                 xlim([450 480])
%                 
%                 figure;
%                 tiledlayout(1,2)
%                 nexttile
%                 plot(bgCh3Temp)
%                 xlim([450 480])
%                 nexttile
%                 plot(obj.bgCh3Aligned)
%                 xlim([450 480])
                %-----------------------------------------------average data--------------------------------------------------------------------------
                bg1 = mean(obj.bgCh1Aligned,2);
                bg2 = mean(obj.bgCh2Aligned,2);
                bg3 = mean(obj.bgCh3Aligned,2);
                obj.bgCh1 = bg1;
                obj.bgCh2 = bg2;
                obj.bgCh3 = bg3;
                obj.CtrlV1 = output.Channel_1.CtrlV.data(1); % get APD control voltage, use only first one as all are the same
                obj.CtrlV2 = output.Channel_2.CtrlV.data(1); % get APD control voltage, use only first one as all are the same
                obj.CtrlV3 = output.Channel_3.CtrlV.data(1); % get APD control voltage, use only first one as all are the same
            catch
                bg1 = output.Background.Ch1.data;
                bg2 = output.Background.Ch2.data;
                bg3 = output.Background.Ch3.data;
                obj.bgCh1 = bg1;
                obj.bgCh2 = bg2;
                obj.bgCh3 = bg3;
                obj.CtrlV1 = NaN;
                obj.CtrlV2 = NaN;
                obj.CtrlV3 = NaN;
            end
        end
    end
    
end