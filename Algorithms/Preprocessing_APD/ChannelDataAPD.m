classdef ChannelDataAPD < handle
    % top level class for deconvolution (Laguerre, multi-exponential and phasor) of FLImBRUSH data
    % See also APDCLASS, FLIMDATACLASS, SYSCALIDATACLASS.

    properties
        softwareVersion %
        averagedData % averaged data
        alpha % alpha value
        APDObj % apd detector class obj, contain gain and iRF information
        badDataIdx % raw waveform bad data index, result of saturation filtering
        bg % channel background, DC removed
        bgAligned % channel background, DC removed, aligned with data
        bgLow; % BG removal index low
        bgHigh; % BG removal index low
        bg_Estimated; % estimated BG from data
        bw = 400 % MHz detector bandwidth
        CtrlV % APD control voltage
        BgEstimated = 0; % flag for whether ch1 bg is estimated from data
        dtRaw % raw data dt
        dtUp % upsampled data dt
        dataAveragingLV % maximum number of WF to average
        dataAveraging % data averaging used for processing
        dataLow; % data amplitude threshold low
        dataHigh; % data amplitude threshold high
        %         dcLow; % DC removal index low
        %         dcHigh; % DC removal index low
        dataT % truncated data
        goodDataIdx % raw waveform good data index, result of saturation filtering
        gain % gain lise 1D vector
        irfIdx % 1 by (number of data points) index of correponding irf in the APD Obj
        irfGainMax % max gain for variable irf
        irfVMax % max V for variable irf 
        Lg_INTs % intensities
        Lg_INTsGainCorrected % intensities
        K % Laguerre order
        LaguerreBasis = []; % Laguerre base funciton
        Lg_LCs %Laguerre coefficient
        Lg_LTs %lifetimes
        laserRepRate % laser rep rate
        M % truncated data length
        numOfWFs % number of waveforms
        numOfAvgWFs % number of waforms after averaging
        noise % data noise
        outlierFlag % 1D column vector
        preProcessedData % final preprocessed data
        prePeakPoints = 68; % data points before peak
        rawData % raw data
        rawDataDCRemoved % raw data DC removed
        rawDataUpsampled % upsampled raw data after DC removal
        rawCtrlV % raw control voltage, same as number of waveforms
        rawGain % raw gain, same size as number of waveforms
        SNR % data SNR
        stat_test % statistic test
        shift % WF shift amount
        shiftIteration % iteration of recursive search
        Lg_exitflag % flag for convergence, 1: function converge; 0: max iteration achieved
        % shift_fit % fitted shift
        truncationLength % data truncation length
        %         timeStamp % time stamp of averaged data used for image reconstruction
        %         timeStampDecon % time stamp of deconvolved averaged data used for image reconstruction
        upSampleFactor; % upsample factor
        wfLenght % length of each waveform
        wf_aligned % aligned waveform
        exclude % A vector of integers indexing the points you want to exclude, e.g., [1 10 25].
        exp1DeconObj  % single exponential decon object, does not include data that is filtered
        exp2DeconObj  % single exponential decon object, does not include data that is filtered
        exp3DeconObj  % single exponential decon object, does not include data that is filtered
        Ph_H1S % phasor result harmonics 1
        Ph_H1G % phasor result harmonics 1
        Ph_H2S % phasor result harmonics 2
        Ph_H2G % phasor result harmonics 2
        Ph_H3S % phasor result harmonics 3
        Ph_H3G % phasor result harmonics 3
        Ph_H4S % phasor result harmonics 4
        Ph_H4G % phasor result harmonics 4
        %         mExp_a1 % multi-exponential fit result
        %         mExp_a2 % multi-exponential fit result
        %         mExp_a3 % multi-exponential fit result
        %         mExp_tau1 % multi-exponential fit result
        %         mExp_tau2 % multi-exponential fit result
        %         mExp_tau3 % multi-exponential fit result
        %         mExp_LT_dacay % multi-exponential fit result
        %         mExp_LT_formula % multi-exponential fit result
        %         mExp_INT % multi-exponential fit result
    end
    methods
        function obj = ChannelDataAPD(rawDataIn, CtrlVIn, LVAvgIn, APDObjIn, dtIn, bgIn, laserRepRateIn, upSampleFactorIn) % constructor

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
            obj.rawGain = interp1(APDObjIn.gainV,APDObjIn.apdGain,CtrlVIn,'pchip'); % use spline interpolation method
            obj.numOfWFs = size(rawDataIn,2);
            obj.wfLenght = size(rawDataIn,1);
            obj.laserRepRate = laserRepRateIn;
            obj.upSampleFactor = upSampleFactorIn; %updated 09/09/2022
            % calculate SNR
        end

        function getSoftwareVersion(obj) % function to read software version
            ini = IniConfig();
            ini.ReadFile('SoftwareVersion.ini');
            sections = ini.GetSections();
            [keys, ~] = ini.GetKeys(sections{1});
            temp = ini.GetValues(sections{1}, keys);
            MAJOR = temp{1};
            MINOR = temp{2};
            PATCH = temp{3};
            obj.softwareVersion = [num2str(MAJOR) '.' num2str(MINOR) '.' num2str(PATCH)];
        end

        function removeDCData(obj,varargin)
            % switch nargin
            %     case 1
            %         WF_max_threshold = 1.6;
            %     case 2
            %         WF_max_threshold = varargin{1};
            % end
            % numOfWF = size(obj.rawData,2);
            % DC = zeros(1,numOfWF);
            % WFWindow = 500; % WF window to average DC
            % gainWindow = 5; % gain window to average DC
            % satIndex = find(max(obj.rawData)>WF_max_threshold); % saturated waveform index,
            % for i = 1:numOfWF %loop thorugh all data and find DC
            %     G = obj.rawGain(i);
            %     timeIdx1 = i-WFWindow;
            %     if timeIdx1<1
            %         timeIdx1 = 1;
            %     end
            %     timeIdx2 = i+WFWindow;
            %     if timeIdx2 > numOfWF
            %         timeIdx2 = numOfWF;
            %     end
            %     timeIdx =  timeIdx1:timeIdx2; % idx of data point inside time window
            %     Acommon = intersect(timeIdx,satIndex); % find saturated index if it is in time window
            %     timeIdx = setxor(timeIdx,Acommon); % remove sturated waveform index
            %     GTemp = obj.rawGain(timeIdx); % gain of data inside time window
            %     timeIdx(~(GTemp>=G-gainWindow&GTemp<=G+gainWindow)) = []; % remove index of gain out of range
            %     DCAllWF = obj.rawData(end-50:end,timeIdx); % get all data used for DC removal from tail
            %     % DCAllWF = obj.rawData(1:50,timeIdx); % get all data used for DC removal from front
            %     DC(i) = mean(DCAllWF(:)); % average all DC data
            % end
            DC = mean(obj.rawData(end-50:end,:));
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
            obj.rawDataUpsampled = WFUpsampled;
            obj.preProcessedData = WFUpsampled;
            %             figure; plot(WFUpsampled)
            %--------------------------------------------------------------
        end

        function estimateBG(obj)
            %-----------------------estimate BG from data------------------------------
            tempMax = max(obj.rawDataUpsampled);
            tempIdx = find(tempMax>0.8&tempMax<1.1); % find waveform with peak V between 0.8 and 1.1
            if isempty(tempIdx) % if no waveform with max between 0.8 and 1.1
                [~,I] = sort(tempMax, 'descend');
                estimatedBG = obj.rawDataUpsampled(:,I(1:100)); % use highest 100
            else % if there are waveform with max between 0.8 and 1.1
                temp = obj.rawDataUpsampled(:,tempIdx);
                tempCtrlV = obj.rawCtrlV(tempIdx);
                [~,I] = sort(tempCtrlV,'descend'); % sort by gain
                estimatedBG = temp(:,I(1:100)); % use highest 100 gain
            end
            estimatedBGAvg = mean(estimatedBG,2);
            obj.bg_Estimated = estimatedBGAvg;
        end

        function alignWF_CFD(obj, f, range_in) % constant factor alginment of waveforms, each channel has a different factor
            % f is factor
            switch nargin
                case 2 % if no range input, use whole length of data
                    range = 1:size(obj.preProcessedData,1);
                case 3 % if range is specified, use range
                    range = range_in;
                otherwise
                    error('Wrong number of Inputs')
            end
            obj.preProcessedData = alignWaveform_CFDNew(obj.preProcessedData, 2.4, obj.dtUp,f,range);
            %             obj.preProcessedData = alignWaveform_CFDNew(obj.preProcessedData, 2.8, obj.dtUp,f,680:size(obj.preProcessedData,1));
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
            %             nanCol = isnan(sum(obj.preProcessedData));
            %             obj.nonDeconIdx = find(nanCol);
            %             obj.deconIdx = find(~nanCol);
            %             obj.CtrlVDecon = obj.CtrlV(obj.deconIdx);
            %             obj.gainDecon = obj.gain(obj.deconIdx);
            %             sum(isnan(sum(obj.preProcessedData)))
            %             obj.preProcessedData(:,nanCol)=[]; % remove non-decon data
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

        function removeDCBG(obj, bgLowIn, bgHighIn, channel, varargin)
            switch nargin
                case 4
                    threshold = 0.01;
                case 5
                    threshold = varargin{1};
            end
            if channel == 1
                range = bgLowIn:bgHighIn;
            else
                range = bgHighIn:bgHighIn+680-1;
            end
            obj.bgLow = bgLowIn;
            obj.bgHigh = bgHighIn;
            
            dataBGAUC = sum(obj.preProcessedData(bgLowIn:bgHighIn,:));
            idxTemp = find(dataBGAUC>threshold*(bgHighIn-bgLowIn+1)); % WF with BG platau
            dataWithBG = obj.preProcessedData(:,idxTemp);
            if (channel==1)&&(sum(obj.bg_Estimated(range))>1.1*sum(obj.bg(range))) % if channel 1 and estimated BG has higher plateau
                obj.BgEstimated = 1; % use estimated bg
                bgTemp = obj.bg_Estimated;
                warning("Measured BG likely wrong, use estimated BG instead.")
                f = warndlg("Measured BG likely wrong, use estimated BG instead.","BG warning");
                pause(5)
                try
                close(f)
                catch
                end
            else
                obj.BgEstimated = 0;
                bgTemp = obj.bg;
            end
            dataWithBG = [dataWithBG bgTemp]; % attach bg as last waveform
            dataWithBG = alignWaveform_CFDNew(dataWithBG, 2.4, obj.dtUp,0.5,range); % align bg with data
            obj.bgAligned = dataWithBG(:,end); % save aligned bg
            dataWithBG(:,end) = [];
            obj.preProcessedData(:,idxTemp) = dataWithBG;

            dataBGAUC = sum(obj.preProcessedData(bgLowIn:bgHighIn,:));
            dataBGAUC(dataBGAUC<threshold*(bgHighIn-bgLowIn+1)) = 0; % if average BG less than 10 mV, do not do BG subtraction, set factor to 0
            bgAUC = sum(obj.bgAligned(bgLowIn:bgHighIn));
            bgScaleFactor = dataBGAUC./bgAUC;
            %             bgScaleFactor(bgScaleFactor<5) = 0;
            bgScaleFactor(isnan(bgScaleFactor)) = 0;
            obj.preProcessedData = obj.preProcessedData-obj.bgAligned*bgScaleFactor;
            obj.preProcessedData = alignWaveform_CFDNew(obj.preProcessedData, 2.4, obj.dtUp,0.7,bgHighIn:bgHighIn+680-1); % realign BG subtracted WF
            obj.noise = std(obj.preProcessedData(1:500,:),1); % use last first 500 points since data is upsampled
            WFMax = max(obj.preProcessedData);
            obj.SNR = 20*log10(WFMax./obj.noise)'; % covert to column vector
        end


        function truncateData(obj, timeWindowIn)
            obj.prePeakPoints = 68;
            %----------------------------------------resampleing irf to match data--------------------------------------------------------------------
            if (obj.APDObj.irfUpSampleddt~=obj.dtUp) % if dt not match, resample irf
                %                 downFactor = obj.dtUp/obj.APDObj.irfUpSampleddt;
                %                 irfL = size(obj.APDObj.irfUpSampled,1);
                %                 irfTemp = zeros(irfL/downFactor,size(obj.APDObj.irfUpSampled,2)); % initialize irf with correct size
                %                 for i = 1:size(obj.APDObj.irfUpSampled,2)
                %                     [~,idx] = max(obj.APDObj.irfUpSampled(:,i));
                %                     xIdx1 = [idx:-downFactor:1];
                %                     xIdx2 = [idx+downFactor:downFactor:irfL];
                %                     xIdx = sort([xIdx1 xIdx2],'ascend');
                %                     irfTemp(:,i) = obj.APDObj.irfUpSampled(xIdx,i);
                %                 end
                tIRF = (0:size(obj.APDObj.irfUpSampled,1)-1)*obj.APDObj.irfUpSampleddt; %time vector for irf
                tData = (0:size(obj.preProcessedData,1)-1)*obj.dtUp; %time vector for data
                irfTemp = interp1(tIRF,obj.APDObj.irfUpSampled,tData); %interpolate irf
                obj.APDObj.irfDecon = irfTemp;
                obj.APDObj.irfdt = obj.dtUp;
            else % if dt match, do not resample use upsampled
                obj.APDObj.irfDecon = obj.APDObj.irfUpSampled;
                obj.APDObj.irfdt = obj.APDObj.irfUpSampleddt;
            end

            obj.preProcessedData(:,find(min(obj.preProcessedData)<-0.05)) = NaN; % remove negative waveform from computation
            [~,dataMaxIAll] = max(obj.preProcessedData);
            dataMaxIAll(dataMaxIAll==1)=NaN;
            dataMaxI = mode(dataMaxIAll);

            % ------------------------peak alignment-----------------------
            % MaxIDiff = dataMaxIAll-dataMaxI;
            % for i = 1:size(obj.preProcessedData,2)
            %     obj.preProcessedData(:,i) = circshift(obj.preProcessedData(:,i),-MaxIDiff(i));
            % end
            %--------------------------------------------------------------

            dataLength = round(timeWindowIn/obj.dtUp);
            obj.truncationLength = dataLength;
            %--------------------------------------------truncate data--------------------------------------------------------------------
            start_idx = dataMaxI-obj.prePeakPoints;
            if start_idx<1
                start_idx=1;
            end
            endIdx = start_idx+dataLength-1;
            if endIdx < size(obj.preProcessedData,1) % if enough data points
                obj.dataT = obj.preProcessedData(start_idx:endIdx,:);
            else % if not enough data points, add zero to the end
                warning('Not enough data to truncate, append with 0!')
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
            % iRFLength = 500; %
            [~,irfMaxI] = max(obj.APDObj.irfDecon); % find irf max
            irfMaxI = mode(irfMaxI);
            %             plot(circshift(obj.APDObj.irfDecon,irfMaxI-mode(irfMaxI),1))
            irf_start_idx = irfMaxI-obj.prePeakPoints;
            %             irf_start_idx = irfMaxI-50;
            if irf_start_idx+iRFLength-1 < size(obj.APDObj.irfDecon,1) % if enough data points
                obj.APDObj.irfTNorm = obj.APDObj.irfDecon(irf_start_idx:irf_start_idx+iRFLength-1,:);
            else
                obj.APDObj.irfTNorm = obj.APDObj.irfDecon(irf_start_idx:end,:); % if not use all data points
            end
            %-------------------------------------set front of rising edge to 0-----------------------------------------------------
            %             obj.APDObj.irfTNorm(1:380,:) = zeros(size(obj.APDObj.irfTNorm(1:380,:)));
            %------------------------------------set tail of irf to be zero--------------------------------------------------------
            irfRealLength = iRFLength; % real irf length
            obj.APDObj.irfTNorm(irfRealLength:end,:) = zeros(size(obj.APDObj.irfTNorm(irfRealLength:end,:)));
            %----------------------------------- pre shift irf so the max align with data----------------------------------------------------------
            %             [~,irfTMaxIdx] = max(obj.APDObj.irfTNorm); % find irf max
            %             for i = 1:length(irfTMaxIdx)
            %                 obj.APDObj.irfTNorm(:,i) = circshift(obj.APDObj.irfTNorm(:,i),-dataMaxIT+irfTMaxIdx(i));
            %             end
            obj.APDObj.irfTNorm = obj.APDObj.irfTNorm./sum(obj.APDObj.irfTNorm); %normalize by AUC
            %             obj.APDObj.irfTNorm = circshift(obj.APDObj.irfTNorm,dataLength*obj.prePeakFactor-50);
            %             obj.APDObj.irfTNorm = obj.APDObj.irfTNorm./max(obj.APDObj.irfTNorm);  % normalize by peak value
            %try CFD alignment
            %             obj.dataT = alignWaveform_CFDNew(obj.dataT,0.26, obj.dtUp);
            %             figure;
            %             hold on
            %             plot(obj.dataT(:,1:100:end))
            %             plot(obj.APDObj.irfT(:,1:10:end))
            %             hold off
        end

        function runDeconLG(obj, exclude_in, varargin) % Laguerre Deconvolution
            obj.getSoftwareVersion();
            %---------------generate laguerre functions------------------------
            obj.dataT = double(obj.dataT);
            obj.wf_aligned = zeros(size((obj.dataT)));
            obj.M = size(obj.dataT,1);
            obj.exclude = exclude_in;
            switch nargin
                case 2 %runDecon(obj, exclude)
                    obj.K = 12; %default Laguerre order
                    obj.alpha = alpha_up(obj.M,obj.K);
                    gainCapIn = 10;
                case 3 % runDecon(obj, exclude, order)
                    obj.K = varargin{1};
                    obj.alpha = alpha_up(obj.M,obj.K);
                    gainCapIn = 10;
                case 4 % runDecon(obj, exclude, order, alpha)
                    obj.K = varargin{1};
                    obj.alpha = varargin{2};
                    gainCapIn = 10;
                case 5 % runDecon(obj, exclude, order, alpha, gainCapIn)
                    obj.K = varargin{1};
                    obj.alpha = varargin{2};
                    gainCapIn = varargin{3};
                otherwise
                    warning('Too many input argument for LaguerreModel constructor!')
            end
            obj.LaguerreBasis = Laguerre(obj.M,obj.K,obj.alpha);
            %------------------Initilize result matrix---------------------
            numOfDataPoints = obj.numOfAvgWFs;
            obj.Lg_LCs = zeros(obj.K,numOfDataPoints);
            obj.Lg_LTs = zeros(numOfDataPoints,1);
            obj.Lg_INTs = zeros(numOfDataPoints,1);
            obj.shift = zeros(numOfDataPoints,1);
            obj.shiftIteration = zeros(numOfDataPoints,1);
            obj.Lg_exitflag = zeros(numOfDataPoints,1);
            %             obj.spec_aligned = zeros(size(obj.dataT));
            %             obj.fit = zeros(size(obj.dataT));
            %             obj.res = zeros(size(obj.dataT));
            obj.irfIdx = ones(numOfDataPoints,1);
            %--------------------------------------------------------------

            irfGain = interp1(obj.APDObj.gainV,obj.APDObj.apdGain,obj.APDObj.irfV,'pchip');
            irfVMaxIdx = find(irfGain<=gainCapIn,1,'last')+1; % find 1st index of irf gain large than gainCapIn
            obj.irfVMax = obj.APDObj.irfV(irfVMaxIdx); % volatge to cap irf
            obj.irfGainMax = irfGain(irfVMaxIdx);
            % numOfiRFV = find(obj.APDObj.irfV>obj.irfVMax,1); % find index of irf V higher than 10
            numOfiRFV = irfVMaxIdx; 
            
            % obj.APDObj.irfTNorm(:,211) = obj.APDObj.irfTNorm(:,208);
            % numOfiRFV = 80;
            % loop through all V and find corresponding data index
            for i = 2:1:numOfiRFV % to account for out of range V
                vLow = obj.APDObj.irfV(i-1);
                if i==numOfiRFV
                    % if i+10>numOfiRFV
                    vHigh = 7;
                else
                    vHigh = obj.APDObj.irfV(i);
                end
                idx = find(obj.CtrlV>=vLow & obj.CtrlV< vHigh);
                if ~isempty(idx)
                    spec = obj.dataT(:,idx); % get waveforms
                    irf = obj.APDObj.irfTNorm(:,i); % get irf
                    obj.irfIdx(idx) = i; %  store irf index
                    channelDataStruct = ChannelData(spec,irf,obj.dtUp,1.5,1:length(idx),std(spec(end-50:end,:),1),obj.gain(idx));
                    laguerreObj = LaguerreModel(channelDataStruct,obj.K, obj.alpha);
                    laguerreObj.warmup_laguerre(obj.exclude);
                    laguerreObj.estimate_laguerre([]);
                    obj.Lg_LCs(:,idx) = laguerreObj.LCs;
                    obj.Lg_LTs(idx) = laguerreObj.LTs;
                    obj.Lg_INTs(idx) = laguerreObj.INTs;
                    obj.shift(idx) = laguerreObj.shift;
                    obj.shiftIteration(idx) = laguerreObj.iteration;
                    obj.Lg_exitflag(idx) = laguerreObj.exitflag;
                    obj.wf_aligned(:,idx) = laguerreObj.WF_aligned;
                    %                     obj.fit(:,idx) = get(laguerreObj,'fit');
                    %                     obj.res(:,idx) = get(laguerreObj,'res');
                    %                     obj.spec_aligned(:,idx) = laguerreObj.spec_aligned;
                end
            end
            obj.stat_test = test_stats(get(obj,'wf_aligned'),get(obj,'fit'), obj.dtUp, 20, {'',{'vars', obj.noise.^2}}); % use bw of 20 so no downsmpling of residual
            obj.Lg_INTsGainCorrected = obj.Lg_INTs./obj.gain;
        end

        function runDeconExp(obj, numOfExp, weight, tauLow, tauHigh, exclude_in) % Multiexponential Deconvolution
            %             wf_aligned = zeros(obj.truncationLength,obj.numOfAvgWFs);
            wf_aligned = get(obj,'wf_aligned'); % truncated and aligned waveform for deconvolution
            %             irf = zeros(size(obj.APDObj.irfTNorm,1),obj.numOfAvgWFs);
            irf = obj.APDObj.irfTNorm(:,obj.irfIdx); % irf matrix
            fieldName = strcat('exp',num2str(numOfExp),'DeconObj');
            fieldName = string(fieldName);
            obj.(fieldName) = ExpModel(numOfExp, wf_aligned, irf, obj.dtUp, weight, tauLow, tauHigh, exclude_in);
            runDecon(obj.(fieldName));
        end

        function runPhasor(obj, Harmonics)
            wf_aligned = get(obj,'wf_aligned');
            irf = obj.APDObj.irfTNorm(:,obj.irfIdx);
            phasorOut = zeros(obj.numOfAvgWFs,1);
            parfor i = 1:obj.numOfAvgWFs
                S = wf_aligned(:,i);
                Ir = irf(:,i);
                phasorOut(i)=ComputePhasor(S,Ir,Harmonics,1);
            end
            switch Harmonics
                case 1
                    obj.Ph_H1G = real(phasorOut);
                    obj.Ph_H1S = imag(phasorOut);
                case 2
                    obj.Ph_H2G = real(phasorOut);
                    obj.Ph_H2S = imag(phasorOut);
                case 3
                    obj.Ph_H3G = real(phasorOut);
                    obj.Ph_H3S = imag(phasorOut);
                case 4
                    obj.Ph_H4G = real(phasorOut);
                    obj.Ph_H4S = imag(phasorOut);
            end
        end

        function toSingle(obj)
            obj.averagedData = single(obj.averagedData);
            obj.preProcessedData = single(obj.preProcessedData);
            obj.rawDataUpsampled = single(obj.rawDataUpsampled);
            obj.rawDataDCRemoved = single(obj.rawDataDCRemoved);
            obj.dataT = single(obj.dataT);
            obj.wf_aligned = single(obj.wf_aligned);
            obj.rawData = single(obj.rawData);
        end

        % method to get properties
        % function to get computed result, like fitting, aligned waveform
        % and residue.
        % the 3rd input can be point idx, if empty outout all points
        function result = get(obj,option,varargin)
            switch option % check option
                %---------------------------------get fitting-----------------------------
                case 'fit'
                    if ~isempty(obj.Lg_LCs) % check laguerre coefficient, if empty, run decon
                        switch nargin % check number of function input
                            case 2
                                result = obj.dataT;
                                for i = 1:obj.numOfAvgWFs
                                    irfTemp = obj.APDObj.irfTNorm(:,obj.irfIdx(i));
                                    result(:,i) = filter(irfTemp,1,obj.LaguerreBasis)*obj.Lg_LCs(:,i);
                                end
                            case 3
                                idx = varargin{1};
                                irfTemp = obj.APDObj.irfTNorm(:,obj.irfIdx(idx));
                                result = filter(irfTemp,1,obj.LaguerreBasis)*obj.Lg_LCs(:,idx);
                        end
                    else
                        warning('Deconvolution result not available, run deconvolution before accessing fitted curve!')
                        result = [];
                    end
                    %-----------------------------aligned waveform-----------------------------
                case 'wf_aligned'
                    if ~isempty(obj.wf_aligned) % check whether wf_aligned is empty
                        switch nargin % check number of function input
                            case 2
                                result = obj.wf_aligned;
                            case 3
                                idx = varargin{1};
                                result = obj.wf_aligned(:,idx);
                        end
                    else
                        warning('Deconvolution result not available, run deconvolution before accessing aligned wavefrom!')
                        result = [];
                    end
                    %---------------------------residue----------------------------------------
                case 'res'
                    if ~isempty(obj.Lg_LCs)
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
                        switch nargin % check number of function input
                            case 2
                                result = obj.dataT;
                                for i = 1:obj.numOfAvgWFs
                                    result(:,i) = circshift(obj.dataT(:,i),obj.shift(i));
                                end
                            case 3
                                idx = varargin{1};
                                result = circshift(obj.dataT(:,idx),obj.shift(idx));
                        end
                    end
                case 'decay'
                    if ~isempty(obj.Lg_LCs)
                        result = obj.LaguerreBasis*obj.Lg_LCs;
                    else
                        warning('use estimate_laguerre before accessing fitted decay!')
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

                case 'GainCorrectedWF'
                    result = obj.dataT./obj.gain';

                otherwise
                    warning('unknown option!')
                    result = [];
            end
        end
    end
end