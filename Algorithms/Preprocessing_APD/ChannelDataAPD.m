classdef ChannelDataAPD < handle
    properties
        averagedData % averaged data
        alpha % alpha value
        APDObj % apd detector class obj, contain gain and iRF information
        badDataIdx % raw waveform bad data index, result of saturation filtering
        bg % channel background
        bgDcRemoved % DC removed bg
        bgLow; % BG removal index low
        bgHigh; % BG removal index low
        bw = 400 % MHz detector bandwidth
        CtrlV % APD control voltage
        CtrlVDecon % decon data ctrl v
        dtRaw % raw data dt
        dtUp % upsampled data dt
        dataAveragingLV % maximum number of WF to average
        dataAveraging % data averaging used for processing
        dataLow; % data amplitude threshold low
        dataHigh; % data amplitude threshold high
        dcLow; % DC removal index low
        dcHigh; % DC removal index low
        deconIdx; % deconvolution index of preprocessed data
        dataT % truncated data
        goodDataIdx % raw waveform good data index, result of saturation filtering
        gain % gain lise 1D vector
        gainDecon % Decon data point gain
        irfIdx % 1 by (number of data points) index of correponding irf in the APD Obj
        INTs % intensities
        INTsAllGainCorrected % all intensity including non-decon data
        INTsGainCorrected % intensities
        K % Laguerre order
        LaguerreBasis = []; % Laguerre base funciton
        LCs %Laguerre coefficient
        LTs %lifetimes
        LTsAll % all lifetims including non-decon data
        laserRepRate % laser rep rate
        M % truncated data length
        numOfWFs % number of waveforms
        numOfAvgWFs % number of waforms after averaging
        nonDeconIdx; % non deconvolution index of processed data
        noise % data noise
        outlierFlag % 1D column vector
        preProcessedData % final preprocessed data
        rawData % raw data
        rawDataDCRemoved % raw data DC removed
        rawDataUpsampled % upsampled raw data after DC removal
        rawCtrlV % raw control voltage, same as number of waveforms
        rawGain % raw gain, same size as number of waveforms
        SNR % data SNR
        stat_test % statistic test
        shift % WF shift amount
        truncationLength % data truncation length
        timeStamp % time stamp of averaged data used for image reconstruction
        timeStampDecon % time stamp of deconvolved averaged data used for image reconstruction
        upSampleFactor; % upsample factor
        wfLenght % length of each waveform
        %         spec_aligned % aligned WF with iRF
        %         fit % fitting result
        %         res % residue
    end
    methods
        function obj = ChannelDataAPD(rawDataIn, CtrlVIn, LVAvgIn, APDObjIn, dtIn, bgIn, laserRepRateIn)
            obj.dtRaw = dtIn;
            %             figure;plot(rawDataIn(:,1:4));
            %             rawDataIn = alignWaveform_CFDNew(rawDataIn,2.8,obj.dtRaw); % CFD raw waveform 1st
            %             rawDataIn = alignWaveform_CFD(rawDataIn,2.8,obj.dtRaw); % old CFD alignment method
            %             figure;plot(rawDataIn(:,1:4));
            % flat blip
%             for i = 1:size(rawDataIn,2)
%             rawDataIn(305:350,i) = ones(350-305+1,1)*mean(rawDataIn(325:350,i));
%             end
            obj.rawData = rawDataIn;
            obj.preProcessedData = rawDataIn;
            obj.rawCtrlV = CtrlVIn;
            obj.dataAveragingLV = LVAvgIn;
            obj.APDObj = APDObjIn;
            obj.dtRaw = dtIn;
            obj.bg = bgIn;
            obj.rawGain = interp1(APDObjIn.gainV,APDObjIn.apdGain,CtrlVIn,'spline'); % use spline interpolation method
            obj.numOfWFs = size(rawDataIn,2);
            obj.wfLenght = size(rawDataIn,1);
            obj.laserRepRate = laserRepRateIn;
            obj.upSampleFactor = 4;
            % calculate SNR
        end
        function removeDCData(obj)
            numOfWF = size(obj.preProcessedData,2);
            DC = zeros(1,numOfWF);
            WFWindow = 500; % WF window to average DC
            gainWindow = 5; % gain window to average DC
            for i = 1:numOfWF %loop thorugh all data and find DC
                G = obj.rawGain(i);
                timeIdx1 = i-WFWindow;
                if timeIdx1<1
                    timeIdx1 = 1;
                end
                timeIdx2 = i+WFWindow;
                if timeIdx2 > numOfWF
                    timeIdx2 = numOfWF;
                end
                timeIdx =  timeIdx1:timeIdx2; % idx of data point inside time window
                GTemp = obj.rawGain(timeIdx); % gain of data inside time window
                timeIdx(~(GTemp>=G-gainWindow&GTemp<=G+gainWindow)) = []; % remove index of gain out of range
                DCAllWF = obj.rawData(end-50:end,timeIdx); % get all data used for DC removal
                DC(i) = mean(DCAllWF(:)); % average all DC data
            end
            obj.rawDataDCRemoved = obj.rawData-DC; % store data 
            obj.preProcessedData = obj.rawData-DC; % store data
            
        end
        function upSampleData(obj)
            obj.dtUp = obj.dtRaw/obj.upSampleFactor; % compute new time resolution
            %-------------------------upsample waveforms-------------------
            WFUpsampled = zeros(obj.wfLenght*obj.upSampleFactor,obj.numOfWFs);
            temp = obj.rawDataDCRemoved;
            upsamplefactor = obj.upSampleFactor;
            parfor i = 1:obj.numOfWFs
                WFUpsampled(:,i) = interp(temp(:,i),upsamplefactor);
            end
%             obj.rawDataUpsampled = WFUpsampled;
            obj.preProcessedData = WFUpsampled;
            %             figure; plot(WFUpsampled)
            %--------------------------------------------------------------
        end
        
        function alignWF_CFD(obj, f) % constant factor alginment of waveforms, each channel has a different factor
            % f is factor
            obj.preProcessedData = alignWaveform_CFDNew(obj.preProcessedData, 2.8, obj.dtUp,f);
%             obj.rawDataUpsampled = alignWaveform_CFDNew(obj.rawDataUpsampled , 2.8, obj.dtUp,f);
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
            obj.CtrlV = averageGianOrV(obj.rawCtrlV,avgIn);
            obj.gain = averageGianOrV(obj.rawGain,avgIn);
            %             obj.preProcessedData = alignWaveform_CFDNew(obj.preProcessedData,3,obj.dtRaw);
            [obj.preProcessedData, obj.outlierFlag] = averageWF(obj.preProcessedData, avgIn); % NaN will presist, will not be removed
            obj.averagedData = obj.preProcessedData;
            % record deconIdx
            nanCol = isnan(sum(obj.preProcessedData));
            obj.nonDeconIdx = find(nanCol);
            obj.deconIdx = find(~nanCol);
            obj.CtrlVDecon = obj.CtrlV(obj.deconIdx);
            obj.gainDecon = obj.gain(obj.deconIdx);
            %             sum(isnan(sum(obj.preProcessedData)))
            obj.preProcessedData(:,nanCol)=[]; % remove non-decon data
            obj.noise = std(obj.averagedData(end-200:end,:),1); % use last 200 points since data is upsampled
            WFMax = max(obj.averagedData);
            obj.SNR = 20*log10(WFMax./obj.noise)'; % covert to column vector
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
            
%             DC = mean(obj.preProcessedData(dcLowIn:dcHighIn,:));
%             obj.preProcessedData = obj.preProcessedData-DC;
            
            bgDC = mean(obj.bg(dcLowIn:dcHighIn)); % compute BG DC using all BG data, BG is 1D vector
            obj.bgDcRemoved = obj.bg - bgDC;
            
            dataBGAUC = sum(obj.preProcessedData(bgLowIn:bgHighIn,:));
            dataBGAUC(dataBGAUC<3) = 0; % if BG area AUC less than 0.05, do not do BG subtraction, set factor to 0
            bgAUC = sum(obj.bgDcRemoved(bgLowIn:bgHighIn));
            
            bgScaleFactor = dataBGAUC./bgAUC;
            obj.preProcessedData = obj.preProcessedData-obj.bgDcRemoved*bgScaleFactor;
        end
        
        
        function truncateData(obj, timeWindowIn)
            %----------------------------------------resampleing irf to match data--------------------------------------------------------------------
            if (obj.APDObj.irfUpSampleddt~=obj.dtUp) % if dt not match, resample irf
                downFactor = obj.dtUp/obj.APDObj.irfUpSampleddt;
                irfL = size(obj.APDObj.irfUpSampled,1);
                irfTemp = zeros(irfL/downFactor,size(obj.APDObj.irfUpSampled,2)); % initialize irf with correct size
                
                for i = 1:size(obj.APDObj.irfUpSampled,2)
                    [~,idx] = max(obj.APDObj.irfUpSampled(:,i));
                    xIdx1 = [idx:-downFactor:1];
                    xIdx2 = [idx+downFactor:downFactor:irfL];
                    xIdx = sort([xIdx1 xIdx2],'ascend');
                    irfTemp(:,i) = obj.APDObj.irfUpSampled(xIdx,i);
                end
                obj.APDObj.irfDecon = irfTemp;
                obj.APDObj.irfdt = obj.dtUp;
            else % if dt match, do not resample use upsampled
                obj.APDObj.irfDecon = obj.APDObj.irfUpSampled;
                obj.APDObj.irfdt = obj.APDObj.irfUpSampleddt;
            end
            [~,dataMaxI] = max(obj.preProcessedData);
            dataMaxI = mode(dataMaxI);
            
            dataLength = round(timeWindowIn/obj.dtUp);
            obj.truncationLength = dataLength;
            %--------------------------------------------truncate data--------------------------------------------------------------------
            start_idx = dataMaxI-round(0.2*dataLength);
            if start_idx<1
                start_idx=1;
            end
            endIdx = start_idx+dataLength-1;
            if endIdx < size(obj.preProcessedData,1) % if enough data points
                obj.dataT = obj.preProcessedData(start_idx:endIdx,:);
            else % if not enough data points, add zero to the end
                availableWFL = length(obj.preProcessedData(start_idx:end,1));
                temp = zeros(dataLength,size(obj.preProcessedData,2));
                temp(1:availableWFL,:) = obj.preProcessedData(start_idx:end,:);
                obj.dataT = temp;
            end
            %------------------------------------------------------find max idx of truncated data-----------------------------------------------------
%             [~,dataMaxIT] = max(obj.dataT);
%             dataMaxIT = mode(dataMaxIT);
            %----------------------------------------------truncate irf---------------------------------------------------------------------
            iRFLength = dataLength; % use short irf 1000 points
            [~,irfMaxI] = max(obj.APDObj.irfDecon); % find irf max
            irfMaxI = mode(irfMaxI);
            irf_start_idx = irfMaxI-round(0.2*iRFLength);
            if irf_start_idx+iRFLength-1 < size(obj.APDObj.irfDecon,1) % if enough data points
                obj.APDObj.irfTNorm = obj.APDObj.irfDecon(irf_start_idx:irf_start_idx+iRFLength-1,:);
            else
                obj.APDObj.irfTNorm = obj.APDObj.irfDecon(irf_start_idx:end,:); % if not use all data points
            end
            %------------------------------------set tail of irf to be zero--------------------------------------------------------
            irfRealLength = 1000; % real irf length
            obj.APDObj.irfTNorm(irfRealLength:end,:) = zeros(size(obj.APDObj.irfTNorm(irfRealLength:end,:)));
            %----------------------------------- pre shift irf so the max align with data----------------------------------------------------------
%             [~,irfTMaxIdx] = max(obj.APDObj.irfTNorm); % find irf max
%             for i = 1:length(irfTMaxIdx)
%                 obj.APDObj.irfTNorm(:,i) = circshift(obj.APDObj.irfTNorm(:,i),-dataMaxIT+irfTMaxIdx(i)); 
%             end
            obj.APDObj.irfTNorm = obj.APDObj.irfTNorm./sum(obj.APDObj.irfTNorm); %normalize by AUC
            %             obj.APDObj.irfTNorm = obj.APDObj.irfTNorm./max(obj.APDObj.irfTNorm);  % normalize by peak value
            %try CFD alignment
            %             obj.dataT = alignWaveform_CFDNew(obj.dataT,0.26, obj.dtUp);
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
                case 2 % runDecon(obj, order)
                    obj.K = varargin{1};
                    obj.alpha = alpha_up(obj.M,obj.K);
                case 3 % runDecon(obj, order, alpha)
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
            for i = 1:numOfiRFV % to account for out of range V
                if i==1
                    vLow = 0;
                else
                    vLow = obj.APDObj.irfV(i-1);
                end
                if i==numOfiRFV
                    vHigh = 5;
                else
                    vHigh = obj.APDObj.irfV(i);
                end
                idx = find(obj.CtrlVDecon>=vLow & obj.CtrlVDecon< vHigh);
                if ~isempty(idx)
                    spec = obj.dataT(:,idx); % get waveforms
                    irf = obj.APDObj.irfTNorm(:,i); % get irf
                    obj.irfIdx(idx) = i; %  store irf index
                    channelDataStruct = ChannelData(spec,irf,obj.dtUp,1.5,1:length(idx),std(spec(end-50:end,:),1),obj.gainDecon(idx));
                    laguerreObj = LaguerreModel(channelDataStruct,obj.K, obj.alpha);
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
            obj.INTsGainCorrected = obj.INTs./obj.gainDecon;
            obj.LTsAll = zeros(size(obj.averagedData,2),1);
            obj.LTsAll(obj.deconIdx) = obj.LTs;
            obj.INTsAllGainCorrected = zeros(size(obj.averagedData,2),1);
            obj.INTsAllGainCorrected(obj.deconIdx) = obj.INTsGainCorrected;
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
                                    irfTemp = obj.APDObj.irfTNorm(:,obj.irfIdx(i));
                                    result(:,i) = filter(irfTemp,1,obj.LaguerreBasis)*obj.LCs(:,i);
                                end
                            case 3
                                idx = varargin{1};
                                irfTemp = obj.APDObj.irfTNorm(:,obj.irfIdx(idx));
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
                    %-------------------------------------upsample data----------------------------------------------------------
                case 'rawDataUpsampled'
                    WFUpsampled = zeros(obj.wfLenght*obj.upSampleFactor,obj.numOfWFs);
                    temp = obj.rawDataDCRemoved;
                    upsamplefactor = obj.upSampleFactor;
                    parfor i = 1:obj.numOfWFs
                        WFUpsampled(:,i) = interp(temp(:,i),upsamplefactor);
                    end
                    result = WFUpsampled;
                
                otherwise
                    warning('unknown option!')
                    result = [];
            end
        end
    end
end