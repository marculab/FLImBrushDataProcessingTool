classdef ExpModel < handle
    properties (Access = protected)
        channeldata
    end
    
    properties (Access = private)
        M % number of points
        N %number of exponential base
        % initial condition, dummy parameter for now
        initial_cond
    end
    
    properties
        Taus
        Weights
        LTs
        INTs
        stat_test
    end
    
    methods
        % constructor
        function obj = ExpModel(channeldata,N,varargin)
            % varargin: initial_cond
            obj.channeldata = channeldata;
            obj.N = N;
            obj.M = size(channeldata.data,1);
            switch nargin
                case 2
                    obj.initial_cond = [];
                case 3
                    obj.initial_cond = varargin{1};
                otherwise
                    warning('Too many input argument for ExpModel constructor!')
            end
            
        end
        % need to use the iIRF_align function form LaguerreModel
        function iIRF_align(obj)
            %% why????? shift = spec_laser_align(obj.channeldata.data(:,2:end),obj.channeldata.iIRF,8,40,40,[],[])
            shift = spec_laser_align(obj.channeldata.data,obj.channeldata.iIRF,8,40,40,[],[]);
            obj.channeldata.iIRF = circshift(obj.channeldata.iIRF,[shift,0]);
        end
        % Exponential fit
        function fit_exp(obj)
            spec = obj.channeldata.data;
            [obj.Weights, obj.Taus, obj.LTs, obj.INTs, fitt, raw] = ...
                multiexp_fit(spec, obj.channeldata.dt, obj.channeldata.iIRF, obj.N);
            obj.stat_test = test_stats(fitt, raw, obj.channeldata.dt, obj.channeldata.bw);
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
                    if ~isempty(obj.Weights)
                        result = filter(obj.channeldata.iIRF,1,obj.get('decay'));
%                         result = bsxfun(@times,result,max(obj.channeldata.data,[],1));
%                         result = result';
                    else
                        warning('use fit_exp before accessing fitted curve!')
                        result = [];
                    end
                case 'decay'
                    if ~isempty(obj.Weights)
                        timebase = (0:obj.M-1)'*obj.channeldata.dt;
%                         exponentials = multiexp_model(obj.Weights, obj.Taus, timebase, sum(obj.channeldata.data));
                        exponentials = multiexp_model(obj.Weights, obj.Taus, timebase, obj.INTs);
                        result = zeros(length(exponentials),obj.M);
                        for i = 1:length(exponentials)
                            result(i,:) = sum(exponentials{i},2);
                        end
                        result = result';
                    else
                        warning('use fit_exp before accessing fitted decay!')
                        result = [];
                    end
                case 'res'
                    if ~isempty(obj.Weights)
                        result = obj.channeldata.data - obj.get('fit');
                    else
                        warning('use fit_exp before accessing fitted decay!')
                        result = [];
                    end
                case 'M'
                    result = obj.M;
                case 'N'
                    result = obj.N;
                case 'exponentials'
                    timebase = (0:obj.M-1)'*obj.channeldata.dt;
                    result = multiexp_model(obj.Weights, obj.Taus, timebase, obj.INTs);
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
                case 'initial_cond'
                    obj.initial_cond = value;
                otherwise
                    warning('unknown option!')
            end
        end
    end
end