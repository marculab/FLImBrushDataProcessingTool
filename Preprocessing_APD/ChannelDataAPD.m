classdef ChannelDataAPD < handle
    properties
        rawData % raw data
        averagedData % averaged data
        preProcessedData % final preprocessed data
        numOfWFs % number of waveforms
        numOfAvgWFs % number of waforms after averaging
        wfLenght % length of each waveform
        goodDataIdx % raw waveform good data index, result of saturation filtering
        badDataIdx % raw waveform bad data index, result of saturation filtering
        CtrlV % APD control voltage
        CtrlVDecon % decon data ctrl v
        gain % gain lise 1D vector
        gainDecon % Decon data point gain
        APDObj % apd detector class obj, contain gain and iRF information
        bg % channel background
        bgDcRemoved % DC removed bg
        dtRaw % raw data dt
        dtUp % upsampled data dt
        bw = 400 % MHz detector bandwidth
        noise % data noise
        SNR % data SNR
        dataAveragingLV % maximum number of WF to average
        dataAveraging % data averaging used for processing
        dataLow; % data amplitude threshold low
        dataHigh; % data amplitude threshold high
        dcLow; % DC removal index low
        dcHigh; % DC removal index low
        bgLow; % BG removal index low
        bgHigh; % BG removal index low
        deconIdx; % deconvolution index of preprocessed data
        nonDeconIdx; % non deconvolution index of processed data
        dataT % truncated data
        K % Laguerre order
        M % truncated data length
        alpha % alpha value
        LaguerreBasis = []; % Laguerre base funciton
        LCs %Laguerre coefficient
        LTs %lifetimes
        LTsAll % all lifetims including non-decon data
        INTs % intensities
        INTsAllGainCorrected % all intensity including non-decon data
        INTsGainCorrected % intensities
        stat_test % statistic test
        shift % WF shift amount
        %         spec_aligned % aligned WF with iRF
        %         fit % fitting result
        %         res % residue
        irfIdx % 1 by (number of data points) index of correponding irf in the APD Obj
        truncationLength % data truncation length
        timeStamp % time stamp of averaged data used for image reconstruction
        timeStampDecon % time stamp of deconvolved averaged data used for image reconstruction
        laserRepRate % laser rep rate
    end
    methods
        function obj = ChannelDataAPD(rawDataIn, CtrlVIn, MaxWFAvgIn, APDObjIn, dtIn, bgIn, laserRepRateIn)
            obj.rawData = rawDataIn;
            obj.preProcessedData = rawDataIn;
            obj.CtrlV = CtrlVIn;
            obj.dataAveragingLV = MaxWFAvgIn;
            obj.APDObj = APDObjIn;
            obj.dtRaw = dtIn;
            obj.bg = bgIn;
            obj.gain = interp1(APDObjIn.gainV,APDObjIn.apdGain,CtrlVIn,'spline'); % use spline interpolation method
            obj.numOfWFs = size(rawDataIn,2);
            obj.wfLenght = size(rawDataIn,1);
            obj.laserRepRate = laserRepRateIn;
            % calculate SNR
        end
        
        function [totalNumOfWf, goodNumOfWf, badNumOfWf]=removeSaturation(obj, low, high)
            obj.dataLow = low;
            obj.dataHigh = high;
            dataMaxV = max(obj.preProcessedData);
            saturationIdx = find(dataMaxV>high);
            lowSignalIdx = find(dataMaxV<low);
            obj.badDataIdx = union(saturationIdx,lowSignalIdx,'sorted');
            obj.goodDataIdx = setdiff(1:size(obj.rawData,2),obj.badDataIdx);
            goodNumOfWf = length(obj.goodDataIdx);
            badNumOfWf = length(obj.badDataIdx);
            totalNumOfWf = goodNumOfWf+badNumOfWf;
            obj.preProcessedData(:,obj.badDataIdx) = NaN;
%             msg = sprintf('Total number of waveforms: %d. \n Waveforms removed: %d. \n Waveforms left: %d.',...
%                 length(obj.badDataIdx)+length(obj.goodDataIdx),length(obj.badDataIdx),length(obj.goodDataIdx));
%             msgbox(msg);
        end
        
        function avgData(obj,avgIn)
            obj.dataAveraging = avgIn;
            [obj.preProcessedData, ~] = averageWF(obj.preProcessedData, avgIn);
            obj.averagedData = obj.preProcessedData;
            obj.noise = std(obj.averagedData(end-50:end,:),1);
            WFMax = max(obj.averagedData);
            obj.SNR = 20*log10(WFMax./obj.noise);
            obj.numOfAvgWFs = obj.numOfWFs/avgIn;
        end
        
        function sortedData = sortData(obj, order)
            dataMaxV = max(obj.preProcessedData);
            [~,I] = sort(dataMaxV,order);
            sortedData = obj.preProcessedData(:,I);
        end
        
        function sortedTruncatedData = sortTruncatedData(obj, order)
            dataMaxV = max(obj.dataT);
            [~,I] = sort(dataMaxV,order);
            sortedTruncatedData = obj.dataT(:,I);
        end
        
        function removeDCBG(obj,dcLowIn, dcHighIn, bgLowIn, bgHighIn)
            obj.dcLow = dcLowIn;
            obj.dcHigh = dcHighIn;
            obj.bgLow = bgLowIn;
            obj.bgHigh = bgHighIn;
            
            DC = mean(obj.preProcessedData(dcLowIn:dcHighIn,:));
            obj.preProcessedData = obj.preProcessedData-DC;
            
            bgDC = mean(obj.bg(dcLowIn:dcHighIn));
            obj.bgDcRemoved = obj.bg - bgDC;
            
            dataBGAUC = sum(obj.preProcessedData(bgLowIn:bgHighIn,:));
            dataBGAUC(dataBGAUC<0.03) = 0;
            bgAUC = sum(obj.bgDcRemoved(bgLowIn:bgHighIn));
            
            bgScaleFactor = dataBGAUC./bgAUC;
            obj.preProcessedData = obj.preProcessedData-obj.bgDcRemoved'*bgScaleFactor;
        end
        
        function upSampleData(obj)
            obj.dtUp = 0.5*obj.dtRaw;
            nanCol = isnan(sum(obj.preProcessedData));
            obj.nonDeconIdx = find(nanCol);
            obj.deconIdx = find(~nanCol);
            WFtoDecon = obj.preProcessedData(:,obj.deconIdx);
            obj.CtrlVDecon = obj.CtrlV(obj.deconIdx);
            obj.gainDecon = obj.gain(obj.deconIdx);
            %-------------------------upsample waveforms-------------------
            WFUpsampled = zeros(length(interp(WFtoDecon(:,1),2)),size(WFtoDecon,2));
            for i = 1:size(WFtoDecon,2)
                WFUpsampled(:,i) = interp(WFtoDecon(:,i),2);
            end
            [~,peakLocation] = max(WFUpsampled);
            peakLocation = mode(peakLocation);
            WFUpsampled = alignWaveform_CFD(WFUpsampled,0.5,peakLocation);
            obj.preProcessedData = WFUpsampled;
            %             figure; plot(WFUpsampled)
            %--------------------------------------------------------------
        end
        
        function truncateData(obj, timeWindowIn)
            [~,dataMaxI] = max(obj.preProcessedData);
            dataMaxI = mode(dataMaxI);
            [~,irfMaxI] = max(obj.APDObj.irfNorm);
            irfMaxI = mode(irfMaxI);
            obj.APDObj.irfNorm = circshift(obj.APDObj.irfNorm,dataMaxI-irfMaxI);
            
            dataLength = timeWindowIn/obj.dtUp;
            obj.truncationLength = dataLength;
            
            start_idx = dataMaxI-round(0.2*dataLength);
            if start_idx<1
                start_idx=1;
            end
            obj.dataT = obj.preProcessedData(start_idx:start_idx+dataLength-1,:);
            obj.APDObj.irfT = obj.APDObj.irfNorm(start_idx:start_idx+dataLength-1,:);
            %             figure;
            %             hold on
            %             plot(obj.dataT(:,1:100:end))
            %             plot(obj.APDObj.irfT(:,1:10:end))
            %             hold off
        end
        
        function runDecon(obj, varargin)
            %---------------generate laguerre functions------------------------
            obj.M = size(obj.dataT,1);
            switch nargin
                case 1
                    obj.K = 12; %default Laguerre order
                    obj.alpha = alpha_up(obj.M,obj.K);
                case 2
                    obj.K = varargin{1};
                    obj.alpha = alpha_up(obj.M,obj.K);
                case 3
                    obj.K = varargin{1};
                    obj.alpha = varargin{2};
                otherwise
                    warning('Too many input argument for LaguerreModel constructor!')
            end
            obj.LaguerreBasis = Laguerre(obj.M,obj.K,obj.alpha);
            %------------------Initilize result matrix---------------------
            numOfDataPoints = length(obj.deconIdx);
            obj.LCs = zeros(obj.K,numOfDataPoints);
            obj.LTs = zeros(numOfDataPoints,1);
            obj.INTs = zeros(numOfDataPoints,1);
            obj.shift = zeros(numOfDataPoints,1);
%             obj.spec_aligned = zeros(size(obj.dataT));
%             obj.fit = zeros(size(obj.dataT));
%             obj.res = zeros(size(obj.dataT));
            obj.irfIdx = zeros(numOfDataPoints,1);
            %--------------------------------------------------------------
            
            numOfiRFV = length(obj.APDObj.irfV);
            
            % loop through all V and find corresponding data index
            for i = 2:numOfiRFV
                vLow = obj.APDObj.irfV(i-1);
                vHigh = obj.APDObj.irfV(i);
                idx = find(obj.CtrlVDecon>=vLow & obj.CtrlVDecon< vHigh);
                if ~isempty(idx)
                    spec = obj.dataT(:,idx); % get waveforms
                    irf = obj.APDObj.irfT(:,i); % get irf
                    obj.irfIdx(idx) = i; %  store irf index
                    channelDataStruct = ChannelData(spec,irf,obj.dtUp,1.5,1:length(idx),std(spec(end-50:end,:),1),obj.gainDecon(idx));
                    laguerreObj = LaguerreModel(channelDataStruct);
                    laguerreObj.estimate_laguerre();
                    obj.LCs(:,idx) = laguerreObj.LCs;
                    obj.LTs(idx) = laguerreObj.LTs;
                    obj.INTs(idx) = laguerreObj.INTs;
                    obj.shift(idx) = laguerreObj.shift;
                    %                     obj.fit(:,idx) = get(laguerreObj,'fit');
                    %                     obj.res(:,idx) = get(laguerreObj,'res');
                    %                     obj.spec_aligned(:,idx) = laguerreObj.spec_aligned;
                end
            end
            %             obj.stat_test = test_stats(obj.spec_aligned,obj.fit, obj.dtUp, obj.bw);
            obj.INTsGainCorrected = obj.INTs./obj.gainDecon';
            obj.LTsAll = zeros(size(obj.averagedData,2),1);
            obj.LTsAll(obj.deconIdx) = obj.LTs';
            obj.INTsAllGainCorrected = zeros(size(obj.averagedData,2),1);
            obj.INTsAllGainCorrected(obj.deconIdx) = obj.INTsGainCorrected';
            obj.LTsAll(obj.LTsAll==0)=NaN;
            obj.INTsAllGainCorrected(obj.INTsAllGainCorrected==0)=NaN;
        end
        
        % method to get properties
        % function to get computed result, like fitting, aligned waveform
        % and residue.
        % the 3rd input can be point idx, if empty outout all points
        function result = get(obj,option,varargin)
            switch option % check option
                %---------------------------------get fitting-----------------------------
                case 'fit'
                    if ~isempty(obj.LCs) % check laguerre coefficient, if empty, run decon
                        switch nargin % check number of function input
                            case 2
                                result = zeros(obj.truncationLength,length(obj.deconIdx));
                                for i = 1:length(obj.deconIdx)
                                    irfTemp = obj.APDObj.irfT(:,obj.irfIdx(i));
                                    result(:,i) = filter(irfTemp,1,obj.LaguerreBasis)*obj.LCs(:,i);
                                end
                            case 3
                                idx = varargin{1};
                                irfTemp = obj.APDObj.irfT(:,obj.irfIdx(idx));
                                result = filter(irfTemp,1,obj.LaguerreBasis)*obj.LCs(:,idx);
                        end
                    else
                        warning('Deconvolution result not available, run deconvolution before accessing fitted curve!')
                        result = [];
                    end
                    %-----------------------------aligned waveform-----------------------------
                case 'wf_aligned'
                    if ~isempty(obj.LCs)
                        switch nargin % check number of function input
                            case 2
                                result = zeros(obj.truncationLength,length(obj.deconIdx));
                                for i = 1:length(obj.deconIdx)
                                    result(:,i) = circshift(obj.dataT(:,i),obj.shift(i));
                                end
                            case 3
                                idx = varargin{1};
                                result = circshift(obj.dataT(:,idx),obj.shift(idx));
                        end
                    else
                        warning('Deconvolution result not available, run deconvolution before accessing fitted curve!')
                        result = [];
                    end
                    %---------------------------residue----------------------------------------
                case 'res'
                    if ~isempty(obj.LCs)
                        switch nargin % check number of function input
                            case 2
                                fitting = get(obj,'fit');
                                dataAligned = get(obj,'wf_aligned');
                                result = dataAligned - fitting;
                            case 3
                                idx = varargin{1};
                                fitting = get(obj,'fit',idx);
                                dataAligned = get(obj,'wf_aligned',idx);
                                result = dataAligned - fitting;
                        end
                    else
                        warning('Deconvolution result not available, run deconvolution before accessing fitted curve!')
                        result = [];
                    end
                otherwise
                    warning('unknown option!')
                    result = [];
            end
        end
    end
end