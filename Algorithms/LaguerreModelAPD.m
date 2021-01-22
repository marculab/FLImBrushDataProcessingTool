classdef LaguerreModelAPD < handle
    % specification for channeldata and iIRF:
    % 1. channeldata was background & DC subtracted, truncated and gain
    % corrected
    % 2. iIRF is unit-integral scaled, truncated to proper length
    % * 10-90 percent before-after peak positioning is suggested for both
    % channeldata and iIRF
    properties (Access = protected)
         %raw data
    end
    
    properties (Access = private)
        LaguerreBasis = []; % Laguerre base funciton
        M %lenght of data
        K = 12; % Lagueere Order
        alpha % Alpha value
    end
    
    properties
        LCs %Laguerre coefficient
        LTs %lifetimes
        INTs % intensities
        stat_test % statistic test
        channeldata % all relevent parameters and data
        shift % laser shift amount
        AUC % area under curve
        threshold % threshold for data 2 by 1
        goodDataIdx % idx of good measurement point
    end
    
    methods
        % constructor
        function obj = LaguerreModel(channeldata,varargin)
            % argument in: channeldata class, alpha
            obj.channeldata = channeldata;
            obj.M = size(channeldata.data,1);
            obj.AUC = sum(obj.channeldata.data,1);
            obj.goodDataIdx = 1:size(channeldata.data,2);
            obj.LTs = zeros(1,size(channeldata.data,2));
            obj.INTs = zeros(1,size(channeldata.data,2));
            
            % use switch if more arguments were needed in future
            switch nargin
                case 1
                    obj.alpha = alpha_up(obj.M,obj.K);
                case 2
                    obj.alpha = varargin{1};
                otherwise
                    warning('Too many input argument for LaguerreModel constructor!')
            end
            obj.LaguerreBasis = Laguerre(obj.M,obj.K,obj.alpha);
        end
        % align iIRF
        function iIRF_align(obj,varargin)
            switch nargin
                case 1
                    obj.shift = spec_laser_align(obj.channeldata.data,obj.channeldata.iIRF,12,20,100,[],[]);  %shift iRF
                case 2
                    obj.shift = varargin{1};
            end
            obj.channeldata.iIRF = circshift(obj.channeldata.iIRF,[obj.shift,0]);
        end
        
        % detect saturation and low signal before deconvolition
        function dataThresholding(obj, thresholdIn)
            obj.threshold = thresholdIn;
%             fullIdx = 1:size(obj.channeldata.data,2);
            badDataIdx = detectSaturationFunction(obj.channeldata.data,thresholdIn(1),thresholdIn(2));
            fullIdx = 1:size(obj.channeldata.data,2);
            obj.goodDataIdx = setdiff(fullIdx,badDataIdx)';
        end
        % do deconvolution
        function estimate_laguerre(obj)
            % data duplicated here, since communication overhead may incur
            % within the parallel for loop if using "obj.channeldata.data".
            % create full result matrix
            
            obj.LCs = zeros(obj.K,size(obj.channeldata.data,2));
   
            spec = obj.channeldata.data(:,obj.goodDataIdx);
            vv=filter(obj.channeldata.iIRF,1,obj.LaguerreBasis);
            D_mat=conv2(eye(size(spec,1)),[1,-3,3,-1],'valid')'*obj.LaguerreBasis;
            % third order forward finite difference derivative  matrix
            % times laguerre basis, accuracy is only 1st order
            D=D_mat;
            H=vv'*vv; %positive definite matrix
            H_chol=chol(inv(H)); %Cholesky decomposition
            C=H_chol*D';
            l1=H_chol*vv';
            lam=zeros(size(D,1),size(spec,2));
            parfor i=1:size(spec,2)
                d=l1*spec(:,i);
                lam(:,i)=lsqnonneg(C,d);
            end
            obj.LCs(:,obj.goodDataIdx)=(vv'*vv)\(vv'*spec-D'*lam);
%             obj.LCs(obj.goodDataIdx,:)=(vv'*vv)\(vv'*spec-D'*lam);
            decays = obj.LaguerreBasis*obj.LCs;
            [obj.LTs,obj.INTs] = h_lifet(decays,obj.channeldata.dt,'average');
            obj.stat_test = test_stats(obj.channeldata.data,obj.get('fit'), obj.channeldata.dt, obj.channeldata.bw);
        end
        
        % get parameters
        function result = get(obj,option)
            switch option
                case 'channeldata'
                    if ~isempty(obj.channeldata)
                        result = obj.channeldata;
                    else
                        warning('No Channel!')
                        result = [];
                    end
                    
                case 'fit'
                    if ~isempty(obj.LCs)
                        result = filter(obj.channeldata.iIRF,1,obj.LaguerreBasis)*obj.LCs;
                    else
                        warning('use estimate_laguerre before accessing fitted curve!')
                        result = [];
                    end
                case 'decay'
                    if ~isempty(obj.LCs)
                        result = obj.LaguerreBasis*obj.LCs;
                    else
                        warning('use estimate_laguerre before accessing fitted decay!')
                        result = [];
                    end
                case 'res'
                    if ~isempty(obj.LCs)
                        result = obj.channeldata.data - obj.get('fit');
                    else
                        warning('use estimate_laguerre before accessing fitted decay!')
                        result = [];
                    end
                case 'M'
                    result = obj.M;
                case 'K'
                    result = obj.K;
                case 'alpha'
                    result = obj.alpha;
                case 'basis'
                    result = obj.LaguerreBasis;
                case 'iRF'
                    result = obj.channeldata.iIRF;
                case 'AUC'
                    result = obj.AUC;
                otherwise
                    warning('unknown option!')
                    result = [];
            end
        end
        % set parameters
        function set(obj,option,value)
            switch option
                case 'K'
                    obj.K = value;
                    obj.LaguerreBasis = Laguerre(obj.M,obj.K,obj.alpha);
                case 'alpha'
                    if isnumeric(value)
                        obj.alpha = value;
                    else
                        obj.alpha = alpha_up(obj.M,obj.K);
                    end
                    obj.LaguerreBasis = Laguerre(obj.M,obj.K,obj.alpha);
                otherwise
                    warning('unknown option!')
            end
        end
    end
end