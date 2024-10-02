classdef LaguerreModel < handle
    % LaguerreModel class, top level object used for Laguerre deconvolution

    properties
        LaguerreBasis = []; % Laguerre base funciton, 2D matrix
        M % lenght of data
        K % Laguerre Order
        alpha % Alpha value of Laguerre base functions
        LCs % Laguerre coefficient, 2D matrix
        LTs % Lifetimes, 1D vector
        INTs % Intensities, 1D vector
        stat_test % statistic test
        channeldataObj % channeldata class object containing raw data and other parameters
        shift % number of data points waveform has to shift to match iRF
        iteration % number of iterations of recursive search
        exitflag % convergence flag
        spec_aligned % aligned waveform
        exclude % index of data excluded from decon due to artifact in the data
        vv % vv matrix with exclude
        vv_full % original vv with no exclude
        D  % D matrix
        l1 % l1 matrix
        C % C matrix
        WF_aligned % waveform aligned with iRF
    end

    methods (Access = public)

        function obj = LaguerreModel(channeldata,varargin)
            % constructor,  create object to run Laguerre deconvolution.
            % Syntax:
            % 1. obj = LaguerreModel(channeldata): create object using Laguerre order = 12, and auto-adjust alpha so that the conditional number of
            % Laguerre base functions are 1.01.
            %
            % 2. obj = LaguerreModel(channeldata, order_in): create object using Laguerre order = order_in, and auto-adjust alpha so that the
            % conditional number of Laguerre base functions are 1.01.
            %
            % 3. obj = LaguerreModel(channeldata, order_in, alpha_in): create object using Laguerre order = order_in, and auto-adjust alpha so that
            % the conditional number of Laguerre base functions are 1.01.
            %
            % See also CHANNELDATA.
            obj.channeldataObj = channeldata;
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

        function warmup_laguerre(obj, exclude_in) % function to generate all decon matrix
            obj.exclude = exclude_in; % exclude data for APD system, [540 640]
            obj.vv=filter(obj.channeldataObj.iIRF,1,obj.LaguerreBasis);
            obj.vv_full=filter(obj.channeldataObj.iIRF,1,obj.LaguerreBasis);
            if ~isempty(exclude_in) % if exlude is not empty
                obj.vv(exclude_in,:) = zeros(size(obj.vv(exclude_in,:))); % ignore data in range 540-640
            end
            obj.D=conv2(eye(size(obj.channeldataObj.data,1)),[1,-3,3,-1],'valid')'*obj.LaguerreBasis;
            H=obj.vv'*obj.vv; %positive definite matrix
            H_chol=chol(inv(H)); %Cholesky decomposition
            obj.C=H_chol*obj.D';
            obj.l1=H_chol*obj.vv';
        end

        function estimate_laguerre(obj,exclude_in) % recursive function to find best shift and decon
            LCs = zeros(obj.K,size(obj.channeldataObj.data,2));
            shift = zeros(size(obj.channeldataObj.data,2),1);
            WF_aligned = zeros(size(obj.channeldataObj.data));
            fit = zeros(size(obj.channeldataObj.data));
            iteration = zeros(size(obj.channeldataObj.data,2),1);
            exitflag = zeros(size(obj.channeldataObj.data,2),1);
            l1 = obj.l1;
            C = obj.C;
            vv = obj.vv;
            D = obj.D;
            WF = obj.channeldataObj.data;
            parfor i = 1:size(obj.channeldataObj.data,2)                
                % [rss, obj.LCs(:,i)] = fitData(obj.channeldataObj.data(:,i), 8.7, obj.l1, obj.C, obj.vv, obj.D);
                if ~isnan(sum(WF(:,i))) % if waveform is not all NaN
                [shift(i), LCs(:,i), WF_aligned(:,i), fit(:,i), iteration(i), exitflag(i)] = findShift(WF(:,i), l1, C, vv, D);
                end
            end
            obj.shift = shift;
            obj.LCs = LCs;
            obj.WF_aligned = WF_aligned;
            obj.iteration = iteration;
            obj.exitflag = exitflag;
            decays = obj.LaguerreBasis*obj.LCs;
            [obj.LTs,obj.INTs] = h_lifet(decays,obj.channeldataObj.dt,'average');
            
        end

        function obj_out = estimate_laguerre_old(obj, exclude_in, varargin)
            % run Laguerre deconvolution to estimate Laguerre coefficients
            % syntax:
            % 1. estimate_laguerre(obj) or obj.estimate_laguerre: run Laguerre deconvolution using all data point and default shift range of -20 to 20
            % points.
            %
            % 2. estimate_laguerre(obj, exclude_in): run Laguerre deconvolution with part of the data excluded and default shift range of -20 to 20
            % points.
            % exclude_in: A vector of integers indexing the points you want to exclude, e.g., [1 10 25].
            %
            % 3. estimate_laguerre(obj, exclude_in, shift_range): run Laguerre deconvolution with part of the data excluded and specified shift range.
            % exclude_in: A vector of integers indexing the points you want to exclude, e.g., [1 10 25].
            % shift_range: A vector of integers specifing the shift range, for example: shift_range= -20:20;

            obj.exclude = exclude_in; % exclude data for APD system, [540 640]
            switch nargin
                case 2
                    shift_range=0:15; %default shift range order
                    %                     shift_range= 0; %default shift range order
                case 3
                    shift_range = varargin{1};

                otherwise
                    warning('Too many input argument for function estimate_laguerre!')
            end
            spec_raw = obj.channeldataObj.data;
            spec = spec_raw;
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

            vv=filter(obj.channeldataObj.iIRF,1,LaguerreBasisS);
            if ~isempty(obj.exclude) % if exlude is not empty
                vv(obj.exclude,:) = zeros(size(vv(obj.exclude,:))); % ignore data in range 540-640
            end
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
            %             exitflag = zeros(size(spec,2),1);
            %             options = optimset('TolX',eps);
            parfor i=1:size(spec,2)
                d=l1*spec(:,i);
                [lam(:,i),~,~,myflag(i),~]=lsqnonneg(C,d);
            end
            obj.LCs=(vv'*vv)\(vv'*spec-D'*lam);
            %             fit_all = obj.get('fit');
            fit_all = filter(obj.channeldataObj.iIRF,1,LaguerreBasisS)*obj.LCs; % get fit without blip
            res = spec-fit_all;
            peakError = zeros(size(res,2),1); % get peak error
            parfor m = 1:size(res,2)
                [~,I]=max(spec(:,m));
                startIdx = max(I-20,1);
                endIdx = min(I+5,length(spec(:,m)));
                peakError(m) = sum(res(startIdx:endIdx,m));
            end
            % res = res;
            %             ind1 = sub2ind([length(shift_range),size(spec_raw,2)],14,25);
            %             ind2 = sub2ind([length(shift_range),size(spec_raw,2)],15,25);
            %             figure;plot(res(:,[ind1,ind2]))
            %             figure;plot(spec(:,[ind1,ind2]))
            %             figure;plot(spec(:,ind1),'.-');hold on;plot(fit_all(:,ind1));hold off
            %             figure;plot(spec(:,ind2),'.-');hold on;plot(fit_all(:,ind2));hold off
            %             res = res(40:50,:);
            res_norm = vecnorm(res,2);
            res_norm = reshape(res_norm,[],size(spec_raw,2));
            peakError = reshape(peakError,[],size(spec_raw,2));
            if size(res_norm,1)==1
                best_fit_idx = 1:length(res_norm);
            else
                [~,res_norm_min_idx] = min(res_norm);
                for m = 1:numel(res_norm_min_idx)
                    if (res_norm_min_idx(m)==1)
                        peakErrorTemp = peakError(res_norm_min_idx(m):res_norm_min_idx(m)+1,m);
                        [~,idxTemp] = min(abs(peakErrorTemp));
                        res_norm_min_idx(m) = res_norm_min_idx(m)+idxTemp-1;
                    elseif (res_norm_min_idx(m)==numel(shift_i))
                        peakErrorTemp = peakError(res_norm_min_idx(m)-1:res_norm_min_idx(m),m);
                        [~,idxTemp] = min(abs(peakErrorTemp));
                        res_norm_min_idx(m) = res_norm_min_idx(m)+idxTemp-2;
                    else
                        peakErrorTemp = peakError(res_norm_min_idx(m)-1:res_norm_min_idx(m)+1,m);
                        [~,idxTemp] = min(abs(peakErrorTemp));
                        res_norm_min_idx(m) = res_norm_min_idx(m)+idxTemp-2;
                    end
                    best_fit_idx = res_norm_min_idx+size(res_norm,1)*(1:size(res_norm_min_idx,2))-size(res_norm,1);
                end
            end
            % shiftTemp = shift_v(best_fit_idx); % use one single shift
            % shiftMode = mode(shiftTemp);
            % best_fit_idx = find(shift_v==shiftMode);
            obj.shift = shift_v(best_fit_idx);
            obj.LCs=obj.LCs(:,best_fit_idx);
            obj.spec_aligned = spec(:,best_fit_idx);
            %             fit = obj.get('fit');
            %             idx = 3;
            %             figure;plot(fit(:,idx));hold on;plot(obj.spec_aligned(:,idx),'*-');hold off
            %             figure;plot(obj.spec_aligned(:,1),'.-');hold on;plot(obj.spec_aligned(:,3),'+-');hold off
            %             figure;plot(spec_raw(:,1),'.-');hold on;plot(spec_raw(:,3),'+-');hold off
            decays = obj.LaguerreBasis*obj.LCs;
            % decays(decays<0)=0; % added to fix negative decay
            %             figure;plot(decays);
            %             figure;plot(decays(:,[1,3]))
            %             figure;plot(decays./max(decays));
            %             fit = filter(obj.channeldataObj.iIRF,1,LaguerreBasisS)*obj.LCs;
            %             figure;plot(spec(:,600));hold on;plot(fit(:,600))
            [obj.LTs,obj.INTs] = h_lifet(decays,obj.channeldataObj.dt,'average');
            %             LTe = h_lifet(decays,obj.channeldataObj.dt,'1/e');
            %             obj.stat_test = test_stats(obj.spec_aligned,obj.get('fit'), obj.channeldataObj.dt, obj.channeldataObj.bw);
            obj_out = obj;
        end


        function result = get(obj,option)
            % GET function to retrive object properties
            %
            % Syntax:
            % out = GET(obj, option): get specified properties from object.
            %
            % option list:
            % 'channeldata': retrive channeldataObj property
            % 'fit': get fitting, same size as data
            % 'decay': get fitted decay
            % 'res': get fitting residual
            % 'iRF': get instrument responsed function

            switch option
                case 'channeldata'
                    if ~isempty(obj.channeldataObj)
                        result = obj.channeldataObj;
                    else
                        warning('No Channel!')
                        result = [];
                    end

                case 'fit'
                    if ~isempty(obj.LCs)
                        result = filter(obj.channeldataObj.iIRF,1,obj.LaguerreBasis)*obj.LCs;
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
                    result = obj.channeldataObj.iIRF;
                otherwise
                    warning('unknown option!')
                    result = [];
            end
        end

        function obj_out = set(obj,option,value)
            % functions to set non-public parameters, depreciated
            switch option
                case 'K'
                    obj.K = value;
                    obj.LaguerreBasis = Laguerre(obj.M,obj.K,obj.alpha);
                    obj_out = obj;
                case 'alpha'
                    if isnumeric(value)
                        obj.alpha = value;
                    else
                        obj.alpha = alpha_up(obj.M,obj.K);
                    end
                    obj.LaguerreBasis = Laguerre(obj.M,obj.K,obj.alpha);
                    obj_out = obj;
                otherwise
                    warning('unknown option!')
            end
        end
    end
end

% some helper functins
function [rss, LCs, WF_aligned, fit] = fitData(WF_in, shift, l1, C, vv, D) % compute Laguerre fitting using pre-computed matrix
    % lam=zeros(size(D,1),1);
    WF = circshift(WF_in, round(shift)); % shift closest integer amount
    dx = shift - round(shift);
    x = 1:length(WF);
    x = x';
    xx = x-dx;
    WF_aligned = interp1(x,WF,xx,"pchip"); % shift decimal amount
    d=l1*WF_aligned;
    lam=lsqnonneg(C,d);
    LCs=(vv'*vv)\(vv'*WF_aligned-D'*lam);
    fit = vv*LCs; % get fit without blip
    res = WF_aligned-fit;
    res(fit==0)=0;
    rss = vecnorm(res,2);
end

function [shift, LCs, WF_aligned, fit, iteration, exitflag] = findShift(WF_in, l1, C, vv, D)
    tol = 0.005;
    s0 = 0;
    s1 = 10;
    s3 = 50;
    R = 0.61803399; % golden ratio
    c = 1-R;
    
    s2 = s1+c*(s3-s1);

    f0 = fitData(WF_in, s0, l1, C, vv, D);
    f1 = fitData(WF_in, s1, l1, C, vv, D);
    f2 = fitData(WF_in, s2, l1, C, vv, D);
    f3 = fitData(WF_in, s3, l1, C, vv, D);
    
    iteration = 0;
    exitflag = 1;
    while (abs(s3-s0) > tol) % recursive bracketing
        iteration = iteration+1;
        if (f2<f1)
            s0 = s1;
            s1 = s2;
            s2 = R*s2+c*s3;

            f1 = f2;
            f2 = fitData(WF_in, s2, l1, C, vv, D);
        else
            s3 = s2;
            s2 = s1;
            s1 = R*s1+c*s0;

            f2 = f1;
            f1 = fitData(WF_in, s1, l1, C, vv, D);
        end
        if iteration>=30 % break if iteration is higher than 30 to avoid infite loop
            exitflag = 0;
            break
        end
    end % end of bracketing

    % sprintf('Iteration %d, tolerence %.5f',iteration, abs(s3-s0)) % print tolerence to command window

    if f1<f2
        shift = s1;
    else
        shift = s2;
    end

    [~, LCs, WF_aligned, fit] = fitData(WF_in, shift, l1, C, vv, D); % re-run decon and get output

end