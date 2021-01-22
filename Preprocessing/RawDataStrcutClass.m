classdef RawDataStrcutClass < handle
    
    properties (Access = protected)
        rawdata %RAWDATA! DO NOT CHANGE!!!!!!2D matrix: total number of points(sum of all lines) x number of time points each data point.
        data %processed data after each step. 2D matrix: total number of points(sum of all lines) x number of time points each data point.
        irf % cell array of iRFs
        bg %1D vector of background. (number of time points each data point by 1)
        gain_value %1D matrix of high voltage gain(total number of points by 1)
        pos_values %1D matrix of class position(total number of points by 1)
        timeStamp %1D vector of timestamp(total number of points by 1)
        spaceStamp % spaceStamp from stage, used to interplate data
        bg_curser_position; %bg subtraction range 1 by 2.
        peak_position; % 1 by 4 vector to store the peak position for each channal
        DC_range;%DC subtraction range 1 by 2.
        line_indx; % line index for reconstruction
        ID; % ID of measurement
    end
    
    properties (Access = public)
        filepath % 1D vector of raw data full file path(total number of files by 1)
        filename % 1D vector of raw data file name(total number of files by 1)
        rawdata_dimention;
        decon_idx % index for deconvolved point
        non_decon_idx % index for non-deconvolved point
        bg_scale_factor; % 1D vector of raw data file name(total number of files by 1)
    end
    methods
        % constructor
        function obj = RawDataStrcutClass(filename,filepath, bg, irf)
            [gain_value,time_stamp,space_stamp,pos_values,line_indx,data_in,ID_in]= LoadRawTRIPLEXData(filename,filepath);
            obj.rawdata = -data_in;
            obj.data = -data_in;
            try
                obj.irf = irf;
            catch
                
            end
            obj.bg = bg;
            obj.gain_value = gain_value;
            obj.pos_values = pos_values;
            obj.timeStamp = time_stamp;
            obj.spaceStamp = space_stamp;
            obj.filepath = filepath;
            obj.filename = filename;
            obj.bg_curser_position = []; % vector
            obj.DC_range = [];% range for DC subtraction
            obj.peak_position = [0 0 0 0];% default peak position is all 0, no channel data
            obj.line_indx = line_indx;
            obj.rawdata_dimention = size(data_in);
            obj.ID = ID_in;
            obj.bg_scale_factor=NaN;
            
        end
    end
    
    methods
        % interface to get parameters
        function out = get(obj, option)
            switch option
                case 'ID'
                    if isempty(obj.ID)
                        warning('No ID is available');
                        out = [];
                    else
                        out = obj.ID;
                    end
                    
                case 'rawdata'
                    if isempty(obj.rawdata)
                        warning('No rawdata is available');
                        out = [];
                    else
                        out = obj.rawdata;
                    end
                    
                case 'data'
                    if isempty(obj.data)
                        warning('No data is available');
                        out = [];
                    else
                        out = obj.data;
                    end
                    
                case 'irf'
                    if isempty(obj.irf)
                        warning('No IRF is available');
                        out = [];
                    else
                        out = obj.irf;
                    end
                    
                case 'bg'
                    if isempty(obj.bg)
                        warning('No background is available');
                        out = [];
                    else
                        out = obj.bg;
                    end
                    
                case 'gain_value'
                    if isempty(obj.gain_value)
                        warning('No gain value is available');
                        out = [];
                    else
                        out = obj.gain_value;
                    end
                    
                case 'pos_values'
                    if isempty(obj.pos_values)
                        warning('No position information is available');
                        out = [];
                    else
                        out = obj.pos_values;
                    end
                    
                case 'timeStamp'
                    if isempty(obj.timeStamp)
                        warning('No timeStamp information is available');
                        out = [];
                    else
                        out = obj.timeStamp;
                    end
                    
                case 'spaceStamp'
                    if isempty(obj.timeStamp)
                        warning('No spaceStamp information is available');
                        out = [];
                    else
                        out = obj.spaceStamp;
                    end
                    
                case 'filepath'
                    if isempty(obj.filepath)
                        warning('No filepath information is available');
                        out = [];
                    else
                        out = obj.filepath;
                    end
                    
                case 'filename'
                    if isempty(obj.filename)
                        warning('No filename information is available');
                        out = [];
                    else
                        out = obj.filename;
                    end
                    
                case 'peak_position'
                    if isempty(obj.peak_position)
                        warning('No bg curser low information is available');
                        out = [];
                    else
                        out = obj.peak_position;
                    end
                    
                case 'line_indx'
                    if isempty(obj.line_indx)
                        warning('No bg curser high information is available');
                        out = [];
                    else
                        out = obj.line_indx;
                    end
                    
                case 'BG_curser'
                    if isempty(obj.bg_curser_position)
                        warning('No bg curser high information is available');
                        out = [];
                    else
                        out = obj.bg_curser_position;
                    end
            end
        end
        
        %call function to remove DC offset in the data
        function out = removeDC(obj,DC_range)
            if isempty(DC_range)
                h = msgbox('Please specify DC subtration range in the menue bar!');
            else
                obj.DC_range = DC_range;
                DC = mean(obj.data(:,DC_range(1):DC_range(2)),2);
                obj.data = bsxfun(@minus, obj.data, DC);
                out = obj.data;
                obj.bg = obj.bg-mean(obj.bg(DC_range(1):DC_range(2)));
            end
        end
        
        % detect saturation bad data is deleted from the matrix
        function [out,non_decon_idx,full_data_size] = detectSaturation(obj,low,high)
            full_data_size = obj.rawdata_dimention;
            non_decon_idx = detectSaturationFunction(obj.data,low,high);
            fullIdx = 1:obj.rawdata_dimention(1);
            fullIdx=fullIdx';
            obj.decon_idx = setdiff(fullIdx,non_decon_idx);
            %             obj.decon_idx = fullIdx;
            %             bad data is deleted from the matrix
            obj.data(non_decon_idx,:)=[];
            % bad data is set to zeros
            %             obj.data(non_decon_idx,:)=zeros(length(non_decon_idx),obj.rawdata_dimention(2));
            out = obj.data;
            obj.non_decon_idx = non_decon_idx;
        end
        % call function to do BG subtraction
        function [out,scale_factor] = removeBG(obj,start_pos,end_pos)
            BG_area = sum(obj.bg(start_pos:end_pos));
            data_area = sum(obj.data(:,start_pos:end_pos),2);
            xx = obj.bg(start_pos:end_pos)';
            yy = obj.data(:,start_pos:end_pos);
            scale_factor_fit=zeros(size(yy,1),1);
            scale_intercept = zeros(size(yy,1),1);
            
            myfittype = fittype('a*x','dependent',{'y'},'independent',{'x'},...
            'coefficients',{'a'});
            xx = conv(xx,ones(1,10)/10,'same');
            yy = conv2(yy,ones(1,10)/10,'same');
%             h10=figure; 
            for i = 1:size(yy,1)
                if mean(yy(i,:))<=0.015
                    scale_factor_fit(i) = 0;
                else
                    idx = yy(i,:)>=-100;
                    P = fit(xx(idx)',yy(i,idx)',myfittype);
                    scale_factor_fit(i) = P.a;
%                     scale_intercept(i) = P(2);
%                     scatter(xx(idx),yy(i,idx))
%                     hold onclc
%                     plot(xx(idx),xx(idx)*P.a)
%                     hold off
%                     xlim([0 max(xx)]);
%                     ylim([0 0.3])
%                     drawnow
%                     pause(0.1)
                end
            end
%             close(h10)
            
            if BG_area==0
                scale_factor = ones(size(obj.data,1),1);
                %                 uiwait(msgbox('BG subrtaction area is 0!, Scale factor is set to 1',...
                %                     'Warning','warning','modal''modal'))
            else
                scale_factor  = data_area/BG_area;
            end
%             h1=figure;plot((obj.data-scale_factor*obj.bg')');
%             pause(5)
%             close(h1)
            scale_factor = scale_factor_fit;
            %             obj.data=obj.data-scale_factor*obj.bg';
            obj.data=obj.data-scale_factor*obj.bg';
            %             ripples = mean(obj.data(1:20,:),1);
            %             obj.data=obj.data-repmat(ripples,length(scale_factor),1);
            obj.bg_curser_position = [start_pos end_pos];
%             obj.bg_scale_factor = scale_factor;
            obj.bg_scale_factor=scale_factor;
            out = obj.data;
            
        end
        
        % function to calculate noise
        function out = get_noise(obj)
            out = std(obj.data(:, (end - 100):end), 0, 2)';
%             out = std(obj.data(:, 1:500), 0, 2)';
        end
        
        % call function to truncate data to 4 channels
        function out = truncateData(obj,truncation_length,peak_position)
            % if peak position=0, no data will be returned for that channel
            peak_left = round(truncation_length*0.1);
            obj.peak_position = peak_position; % record peak position
            % Channel 1
            if ne(peak_position(1), 0)
                
                out{1}=obj.data(:,peak_position(1)-peak_left:peak_position(1)+truncation_length-peak_left-1);
                
            else
                
                out{1} = [];
                
            end
            % Channel 2
            if ne(peak_position(2), 0)
                
                out{2}=obj.data(:,peak_position(2)-peak_left:peak_position(2)+truncation_length-peak_left-1);
                
            else
                
                out{2} = [];
                
            end
            % Channel 3
            if ne(peak_position(3), 0)
                
                out{3}=obj.data(:,peak_position(3)-peak_left:peak_position(3)+truncation_length-peak_left-1);
                
            else
                
                out{3} = [];
                
            end
            % Channel 4
            if ne(peak_position(4), 0)
                
                out{4}=obj.data(:,peak_position(4)-peak_left:peak_position(4)+truncation_length-peak_left-1);
                
            else
                
                out{4} = [];
                
            end
            
            return
            
        end
        
        % call function to do gain correction
        function [gainfactor,gain_matrix] = gainCorrection(obj)
            base_gain = 1500;
            gainfactor = base_gain - obj.gain_value;
            gainfactor = gainfactor / 100;
            gainfactor = 3.^gainfactor;
            base_gain_vector = ones(size(obj.gain_value))*base_gain;
            gain_matrix = [obj.gain_value base_gain_vector];
            
            
        end
        
        
    end
end