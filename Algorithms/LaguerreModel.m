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
        
    end
    
    properties
        LaguerreBasis = []; % Laguerre base funciton
        M %lenght of data
        K % Laguerre Order
        alpha % Alpha value
        LCs %Laguerre coefficient
        LTs %lifetimes
        INTs % intensities
        stat_test % statistic test
        channeldata % all relevent parameters and data
        shift % WF shift amount
        spec_aligned % aligned WF with iRF
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
        %         function iIRF_align(obj,varargin)
        %             switch nargin
        %                 case 1
        %                     n_align = ceil(0.25*size(obj.channeldata.data,2));
        %                     obj.shift = spec_laser_align(obj.channeldata.data,obj.channeldata.iIRF,12,500,10,[],[]);  %shift iRF
        %                 case 2
        %                     obj.shift = varargin{1};
        %             end
        %             obj.channeldata.iIRF = circshift(obj.channeldata.iIRF,[obj.shift,0]);
        %         end
        % do deconvolution
        function estimate_laguerre(obj,varargin)
            % data duplicated here, since communication overhead may incur
            % within the parallel for loop if using "obj.channeldata.data".
            switch nargin
                case 1
                    shift_range= -50:50; %default shift range order
                    %                     shift_range= 0; %default shift range order
                case 2
                    shift_range = varargin{1};
                otherwise
                    warning('Too many input argument for LaguerreModel constructor!')
            end
            spec_raw = obj.channeldata.data;
            spec = spec_raw;
%             spec(585:750,:) = zeros(size(spec(585:750,:))); % remove blip from data
            LaguerreBasisS = obj.LaguerreBasis;
            wfLength = size(spec,1);
            %             spec = spec./max(spec);
            shift_i = shift_range';
            spec = repmat(spec,length(shift_i),1);
            spec = reshape(spec,wfLength,[]);
            shift_v = repmat(shift_i,size(spec_raw,2),1);
            for ii = 1:length(shift_v)
                spec(:,ii) = circshift(spec(:,ii),shift_v(ii));
            end
            
            vv=filter(obj.channeldata.iIRF,1,LaguerreBasisS);
            vv(540:640,:) = zeros(size(vv(540:640,:)));
            D_mat=conv2(eye(size(spec,1)),[1,-3,3,-1],'valid')'*LaguerreBasisS;
%             D_mat(581:585,:) = zeros(size(D_mat(581:585,:)));
            % third order forward finite difference derivative  matrix
            % times laguerre basis, accuracy is only 1st order
            D=D_mat;
            H=vv'*vv; %positive definite matrix
            H_chol=chol(inv(H)); %Cholesky decomposition
            C=H_chol*D';
            l1=H_chol*vv';
            lam=zeros(size(D,1),size(spec,2));
            %             options = optimset('Display','notify','TolX',10*eps);
            exitflag = zeros(size(spec,2),1);
            parfor i=1:size(spec,2)
                d=l1*spec(:,i);
                [lam(:,i),~,~,exitflag(i),~,~]=lsqnonneg(C,d);
            end
            obj.LCs=(vv'*vv)\(vv'*spec-D'*lam);
            %             fit_all = obj.get('fit');
            fit_all = filter(obj.channeldata.iIRF,1,LaguerreBasisS)*obj.LCs; % get fit without blip
            res = spec-fit_all;
            %             ind1 = sub2ind([length(shift_range),size(spec_raw,2)],14,25);
            %             ind2 = sub2ind([length(shift_range),size(spec_raw,2)],15,25);
            %             figure;plot(res(:,[ind1,ind2]))
            %             figure;plot(spec(:,[ind1,ind2]))
            %             figure;plot(spec(:,ind1),'.-');hold on;plot(fit_all(:,ind1));hold off
            %             figure;plot(spec(:,ind2),'.-');hold on;plot(fit_all(:,ind2));hold off
            %             res = res(40:50,:);
            res_norm = vecnorm(res);
            res_norm = reshape(res_norm,[],size(spec_raw,2));
            if size(res_norm,1)==1
                best_fit_idx = 1:length(res_norm);
            else
                [~,res_norm_min_idx] = min(res_norm);
                best_fit_idx = res_norm_min_idx+size(res_norm,1)*(1:size(res_norm_min_idx,2))-size(res_norm,1);
            end
            shiftTemp = shift_v(best_fit_idx); % fixed shift
            shiftMode = mode(shiftTemp);
            best_fit_idx = find(shift_v==shiftMode);
            obj.shift = shift_v(best_fit_idx);
            obj.LCs=obj.LCs(:,best_fit_idx);
            obj.spec_aligned = spec(:,best_fit_idx);
            %             fit = obj.get('fit');
            %             idx = 3;
            %             figure;plot(fit(:,idx));hold on;plot(obj.spec_aligned(:,idx),'*-');hold off
            %             figure;plot(obj.spec_aligned(:,1),'.-');hold on;plot(obj.spec_aligned(:,3),'+-');hold off
            %             figure;plot(spec_raw(:,1),'.-');hold on;plot(spec_raw(:,3),'+-');hold off
            decays = obj.LaguerreBasis*obj.LCs;
            %             figure;plot(decays);
            %             figure;plot(decays(:,[1,3]))
            %             figure;plot(decays./max(decays));
%             fit = filter(obj.channeldata.iIRF,1,LaguerreBasisS)*obj.LCs;
%             figure;plot(spec(:,600));hold on;plot(fit(:,600))
            [obj.LTs,obj.INTs] = h_lifet(decays,obj.channeldata.dt,'average');
            %             LTe = h_lifet(decays,obj.channeldata.dt,'1/e');
            %             obj.stat_test = test_stats(obj.spec_aligned,obj.get('fit'), obj.channeldata.dt, obj.channeldata.bw);
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
                        result = obj.spec_aligned - obj.get('fit');
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