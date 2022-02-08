classdef ExpModel < handle
    % multiexponential deconvolution class
    
    properties
        order = 1; % Laguerre base funciton
        spec_aligned % input waveform aligned with irf
        irf % input instrument response function, column vector or matrix. If matrix, the second dimension must be the same as data
        dt % time resolution
        weight % weight, same length as data, default is equal weight for all data points
        tauUpperLim % upper limit of lifetime component, default value is 0.1 ns
        tauLowerLim % lower limit of lifetime component, default value is 25 ns
        N % number of waveforms for decon
        L % length of waveform
        FO % fit option
        FT % fit type
        A % pre-exponential factors, number of column is the same size as order
        Tau % lifetime components, number of column is the same size as order
        INTs % intensities
        LTs_decay % lifetimes from decay
        LTs_formula % lifetimes from formula
        rSquare % R-squared (coefficient of determination)
        rMSE % Root mean squared error (standard error)
        exitFlag % exit flag, Describes the exit condition of the algorithm. Positive flags indicate convergence, within tolerances. Zero flags indicate that the maximum number of function evaluations or iterations was exceeded. Negative flags indicate that the algorithm did not converge to a solution.
        numOfIteratoin; %iterations
        InitGuess % parameter starting point
        t % time in ns
        exclude % A vector of integers indexing the points you want to exclude, e.g., [1 10 25].
    end
    
    methods
        % constructor
        function obj = ExpModel(order_in, WF_in, irf_in, dt_in, weight_in, tau_lower_lim, tau_upper_lim, exclude_in)
            obj.order = order_in;
            obj.spec_aligned = WF_in;
            obj.N = size(WF_in,2);
            obj.L = size(WF_in,1);
            
            if (size(irf_in,2)==1)||(size(irf_in,2)==obj.N) % check the size of irf
                obj.irf = irf_in;
            else
                error('Size of irf does not match the size of waveforms!')
            end
            
            obj.dt = dt_in;
            obj.exclude = exclude_in;
            
            if ~isempty(weight_in) % if empty use defaut
                obj.weight = weight_in;
            else
                obj.weight = ones(obj.L,1);
            end
            
            if ~isempty(tau_lower_lim) % if empty use defaut
                obj.tauLowerLim = tau_lower_lim;
            else
                obj.tauLowerLim = 0.1;
            end
            
            if ~isempty(tau_upper_lim) % if empty use defaut
                obj.tauUpperLim = tau_upper_lim;
            else
                obj.tauUpperLim = 25;
            end
            
            obj.A = zeros(obj.N,obj.order);
            obj.Tau = zeros(obj.N,obj.order);
            obj.INTs = zeros(obj.N,1);
            
            obj.rSquare = zeros(obj.N,1);
            obj.rMSE = zeros(obj.N,1);
            obj.exitFlag = zeros(obj.N,1);
            obj.numOfIteratoin = zeros(obj.N,1);
            
            t = linspace(0,obj.L-1,obj.L)*obj.dt;
            t =t'; % covert to column vector
            obj.t = t;
            % creat fitting object and fit option
            FMax = 10; % pre-exponential factor max
            switch obj.order
                case 1
                    lower = [0 obj.tauLowerLim];
                    upper = [FMax obj.tauUpperLim];
                    obj.InitGuess = [0.5 4];
                    obj.FO = fitoptions('Method', 'NonlinearLeastSquares','Lower',lower,'Upper',upper,'StartPoint',obj.InitGuess, 'Exclude', obj.exclude);
                    obj.FT = fittype('monoexp_model(x,a,t,L)','problem','L','options',obj.FO);
                case 2
                    lower = [0 0 obj.tauLowerLim obj.tauLowerLim];
                    upper = [FMax FMax obj.tauUpperLim obj.tauUpperLim];
                    obj.InitGuess = [0.5 0.5 4 4];
                    obj.FO = fitoptions('Method', 'NonlinearLeastSquares','Lower',lower,'Upper',upper,'StartPoint',obj.InitGuess, 'Exclude', obj.exclude);
                    obj.FT = fittype('biexp_model(x,a1,a2,t1,t2,L)','problem','L','options',obj.FO);
                case 3
                    lower = [0 0 0 obj.tauLowerLim obj.tauLowerLim obj.tauLowerLim];
                    upper = [FMax FMax FMax obj.tauUpperLim obj.tauUpperLim obj.tauUpperLim];
                    obj.InitGuess = [0.5 0.5 0.5 4 4 4];
                    obj.FO = fitoptions('Method', 'NonlinearLeastSquares','Lower',lower,'Upper',upper,'StartPoint',obj.InitGuess, 'Exclude', obj.exclude);
                    obj.FT = fittype('triexp_model(x,a1,a2,a3,t1,t2,t3,L)','problem','L','options',obj.FO);
                case 4
                    lower = [0 0 0 0 obj.tauLowerLim obj.tauLowerLim obj.tauLowerLim obj.tauLowerLim];
                    upper = [FMax FMax FMax FMax obj.tauUpperLim obj.tauUpperLim obj.tauUpperLim obj.tauUpperLim];
                    obj.InitGuess = [0.5 0.5 0.5 0.5 4 4 4 4];
                    obj.FO = fitoptions('Method', 'NonlinearLeastSquares','Lower',lower,'Upper',upper,'StartPoint',obj.InitGuess, 'Exclude', obj.exclude);
                    obj.FT = fittype('quadriexp_model(x,a1,a2,a3,a4,t1,t2,t3,t4,L)','problem','L','options',obj.FO);
                otherwise
                    error('Multi-exponential fit with %d orders was not implemented!',obj.order);
            end
        end
        
        function runDecon(obj)
            wfALL = obj.spec_aligned; % copy data to avoid parfor communication overhead
            M = max(wfALL); % find maximum of each waveform
            % if no waveform has peak value higher than 0.1, skip decon and return
            if ~any(M>0.1)
                return
            end
            tt =  obj.t;
            AA = obj.A; % local variable of A
            TT = obj.Tau; % local variable of Tau
            O = obj.order; % copy order
            RSquare = obj.rSquare;
            RMSE = obj.rMSE;
            ExitFlag = obj.exitFlag;
            ITERATION = obj.numOfIteratoin;
            wfAvg = mean(wfALL(:,M>0.1),2); % find the average of waveforms with peak value larger than 0.1
            %---------------------find initial guess ------------------------------------
            % find fastest irf which has the max peak value
            irf_max_value = max(obj.irf);
            [~,I] = max(irf_max_value);
            irf_f = obj.irf(:,I);
            [f,~,~] = fit(tt,wfAvg,obj.FT,'problem',irf_f);
            switch obj.order
                case 1
                    obj.InitGuess = [f.a f.t];
                case 2
                    obj.InitGuess = [f.a1 f.a2 f.t1 f.t2];
                case 3
                    obj.InitGuess = [f.a1 f.a2 f.a3 f.t1 f.t2 f.t3];
                case 4
                    obj.InitGuess = [f.a1 f.a2 f.a3 f.a4 f.t1 f.t2 f.t3 f.t4];
                otherwise
                    error('Unknow initial guess')
            end
            obj.FO.StartPoint = obj.InitGuess; % overright initial guess
            obj.FT = setoptions(obj.FT, obj.FO); % overright fit option
            %----------------------run deconvolution---------------------------------------
            FitObj = obj.FT; % copy fit option to avoid calling obj in parfor
            
            if size(obj.irf,2)==1 % if irf is column vector replicate
                irff_all = repmat(obj.irf, 1, obj.N);
            else
                irff_all = obj.irf;
            end
            
            parfor i = 1:obj.N
                wf = wfALL(:,i);
                irff = irff_all(:,i);
                [M,~] = max(wf);
                if (M > 0.1)
                    [ff,gof,out] = fit(tt,wf,FitObj,'problem',irff);
                    switch O
                        case 1
                            AA(i,:) = ff.a;
                            TT(i,:) = ff.t;
                        case 2
                            AA(i,:) = [ff.a1 ff.a2];
                            TT(i,:) = [ff.t1 ff.a2];
                        case 3
                            AA(i,:) = [ff.a1 ff.a2 ff.a3];
                            TT(i,:) = [ff.t1 ff.t2 ff.t3];
                        case 4
                            AA(i,:) = [ff.a1 ff.a2 ff.a3 ff.a4];
                            TT(i,:) = [ff.t1 ff.t2 ff.t3 ff.t4];
                        otherwise
                            
                    end
                    RSquare(i) = gof.rsquare;
                    RMSE(i) = gof.rmse;
                    ExitFlag(i) = out.exitflag;
                    ITERATION(i) = out.iterations;
                end
            end
            obj.rSquare = RSquare;
            obj.rMSE = RMSE;
            obj.exitFlag = ExitFlag;
            obj.numOfIteratoin = ITERATION;
            %-------------sort result by taus---------------------
            [TT, II] = sort(TT,2);
            %             AAA = AA;
            for k =1: obj.N
                AA(k,:) = AA(k,II(k,:));
            end
            %-----------------------------------------------------------
            AA(AA==0)=NaN;
            obj.A = AA;
            TT(TT==0)=NaN;
            obj.Tau = TT;
            decays = get(obj,'decay');
            [obj.LTs_decay, obj.INTs] =h_lifet(decays,obj.dt,'average');
            %------------------------------get formula lifetime------------------------------------------------------------
            NOMINATOR = zeros(obj.N,1);
            deNOMINATOR = zeros(obj.N,1);
            for n = 1:obj.order
                NOMINATOR = NOMINATOR+AA(:,n).*TT(:,n).^2;
                deNOMINATOR = deNOMINATOR+AA(:,n).*TT(:,n);
            end
            obj.LTs_formula = NOMINATOR./deNOMINATOR;
            
            
            
        end
        
        function out = get(obj,option,varargin)
            switch option
                case 'fit'
                    switch nargin % check number of function input
                        case 2 % all fit
                            decay = get(obj,'decay');
                            fit = zeros(obj.L,obj.N);
                            if size(obj.irf,2)==1
                                fit = filter(obj.irf,1,decay);
                            else
                                for i = 1:obj.N
                                    fit(:,i) = filter(obj.irf(:,i),1,decay(:,i));
                                end
                            end
                           
                        case 3 % 1 fit
                            idx = varargin{1};
                            decay = get(obj,'decay',idx);
                            if size(obj.irf,2)==1
                                fit = filter(obj.irf,1,decay);
                            else
                                fit = filter(obj.irf(:,idx),1,decay);
                            end
                    end
                    out = fit;
                case 'decay'
                    switch nargin % check number of function input
                        case 2 % all decays
                            out = zeros(obj.L,obj.N);
                            AA = obj.A;
                            TauT = obj.Tau;
                        case 3 % single decay
                            idx = varargin{1};
                            out = zeros(obj.L,1);
                            AA = obj.A(idx,:);
                            TauT = obj.Tau(idx,:);
                    end
                    
                    tt = obj.t;
                    
                    switch obj.order
                        case 1
                            for i = 1:size(out,2)
                                out(:,i) = AA(i).*exp(-tt./TauT(i));
                            end
                        case 2
                            for i = 1:size(out,2)
                                out(:,i) = AA(i,1).*exp(-tt./TauT(i,1))+AA(i,2).*exp(-tt./TauT(i,2));
                            end
                        case 3
                            for i = 1:size(out,2)
                                out(:,i) = AA(i,1).*exp(-tt/TauT(i,1))+AA(i,2).*exp(-tt./TauT(i,2))+...
                                    AA(i,3).*exp(-tt./TauT(i,3));
                            end
                        case 4
                            for i = 1:size(out,2)
                                out(:,i) = AA(i,1).*exp(-tt./TauT(i,1))+AA(i,2).*exp(-tt./TauT(i,2))+...
                                    AA(i,3).*exp(-tt./TauT(i,3))+AA(i,4).*exp(-tt./TauT(i,4));
                            end
                        otherwise
                            
                    end
            end
        end
    end
end