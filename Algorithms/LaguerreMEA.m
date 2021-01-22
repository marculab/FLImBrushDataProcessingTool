classdef LaguerreMEA < handle
    % specification for channeldata and iIRF:
    % 1. channeldata was background & DC subtracted, truncated and gain
    % corrected
    % 2. iIRF is unit-integral scaled, truncated to proper length
    % * 10-90 percent before-after peak positioning is suggested for both
    % channeldata and iIRF
    properties (Access = protected)
        channeldata
    end
    
    properties (Access = private)
        LaguerreBasis = [];
        M
        K = 12;
        alpha
        expParam
    end
    
    properties
        tau_list
        weight_list
        Exps
        LC_exps
        LCs
        LTs
        INTs
        stat_test
    end
    
    methods
        % constructor
        function obj = LaguerreMEA(channeldata,varargin)
            % argument in: channeldata class, alpha
            obj.channeldata = channeldata;
            obj.M = size(channeldata.data,1);
            % use switch if more arguments were needed in future
            switch nargin
                case 1
                    obj.alpha = alpha_up(obj.M,obj.K);
                    obj.expParam = [0.1 15 0.01];
                case 2
                    obj.alpha = varargin{1};
                    obj.expParam = [0.1 15 0.01];
                case 3
                    obj.alpha = varargin{1};
                    obj.expParam = varargin{2};
                otherwise
                    warning('Too many input argument for LaguerreModel constructor!')
            end
            obj.LaguerreBasis = Laguerre(obj.M,obj.K,obj.alpha);
            obj.tau_list = obj.expParam(1):obj.expParam(3):obj.expParam(2);
            timebase = (0:obj.M-1)'*channeldata.dt;
            obj.Exps = exp(-timebase*obj.tau_list);
            obj.LC_exps = nnls3c(filter(obj.channeldata.iIRF,1,obj.Exps),obj.LaguerreBasis,obj.channeldata.iIRF);
        end
        % align iIRF
        function iIRF_align(obj)
            %% why????? shift = spec_laser_align(obj.channeldata.data(:,2:end),obj.channeldata.iIRF,8,40,40,[],[])
            shift = spec_laser_align(obj.channeldata.data,obj.channeldata.iIRF,8,40,40,[],[]);
            obj.channeldata.iIRF = circshift(obj.channeldata.iIRF,[shift,0]);
        end
        % do deconvolution
        function estimate_laguerre(obj)
            % data duplicated here, since communication overhead may incur
            % within the parallel for loop if using "obj.channeldata.data".
            spec = obj.channeldata.data;
            vv=filter(obj.channeldata.iIRF,1,obj.LaguerreBasis);
            D_mat=conv2(eye(size(spec,1)),[1,-3,3,-1],'valid')'*obj.LaguerreBasis;
            D=D_mat;
            H=vv'*vv;
            H_chol=chol(inv(H));
            C=H_chol*D';
            l1=H_chol*vv';
            
            lam=zeros(size(D,1),size(spec,2));
            
            parfor i=1:size(spec,2)
%             parfor i=1:size(spec,2)
                
                d=l1*spec(:,i);
                lam(:,i)=lsqnonneg(C,d);
                
            end;
            
            obj.LCs=(vv'*vv)\(vv'*spec-D'*lam);
            decays = obj.LaguerreBasis*obj.LCs;
            [obj.LTs,obj.INTs] = h_lifet(decays,obj.channeldata.dt,'average');
            obj.stat_test = test_stats(obj.channeldata.data,obj.get('fit'), obj.channeldata.dt, obj.channeldata.bw);
        end
        
        % analyze multi-exponential components
        function multiExpAnalysis(obj)
            ExpBase = MultiExpModel(tau_list,ones(size(tau_list)),iIRF{channel},time_res,time_res*channel_size(1));
            [basis, Exps, R1, ~] = ExpBase.projLaguerre(80,1e10,alpha_out{channel});
            for i = 1:size(LCs{channel},2)
                %         fit_options = optimset('TolX',1e-12);
                weight_fit0 = lsqnonneg(R1,LCs{channel}(:,i));
                % normalize weight
                weight_fit0 = weight_fit0.*max(basis*R1)';
                
                
                weight_fit = [weight_fit weight_fit0];
                weight_fit2(:,i) = integratesum(tau_list,tau_list2,weight_fit0)'.*weight_calibration';
            end
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