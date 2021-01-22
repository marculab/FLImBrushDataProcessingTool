classdef LaguerreModel < handle
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
        K % Laguerre Order
        alpha % Alpha value
    end
    
    properties
        LCs %Laguerre coefficient
        LTs %lifetimes
        INTs % intensities
        stat_test % statistic test
        channeldata % all relevent parameters and data
        shift % laser shift amount
    end
    
    methods
        % constructor
        function obj = LaguerreModel(channeldata,varargin)
            % argument in: channeldata class, alpha
            obj.channeldata = channeldata;
            obj.M = size(channeldata.data,1);
            % use switch if more arguments were needed in future
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
        end
        % align iIRF
        function iIRF_align(obj,varargin)
            switch nargin
                case 1
                    n_align = ceil(0.25*size(obj.channeldata.data,2));
                    obj.shift = spec_laser_align(obj.channeldata.data,obj.channeldata.iIRF,12,500,10,[],[]);  %shift iRF
                case 2
                    obj.shift = varargin{1};
            end
            obj.channeldata.iIRF = circshift(obj.channeldata.iIRF,[obj.shift,0]);
        end
        % do deconvolution
        function estimate_laguerre(obj)
            % data duplicated here, since communication overhead may incur
            % within the parallel for loop if using "obj.channeldata.data".
            spec = obj.channeldata.data;
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
            obj.LCs=(vv'*vv)\(vv'*spec-D'*lam);
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