function handles = run_Deconvolution(Decon_param, name_cell, folder_cell,Background,iIRF,handles,show_waveform_flag)
%% get necessary parameter
num_of_batches = length(name_cell);


%% get current date and time and generate saving folder name
t=datetime('now');
formatOut = 'mmm-dd-yyyy HH.MM PM';

if Decon_param.Decon_option == 1 %check decon parameetrs for Laguerre
    alpha_max = alpha_up(Decon_param.ch_width,Decon_param.LG_param.order);
    alpha_max = round(alpha_max,3); % round alpha to 3 digit
    if Decon_param.LG_param.alpha > alpha_max
        Decon_param.LG_param.alpha = alpha_max;
        message = sprintf(['Input alpha value out of upper bound, maximum ' ...
            'value of %g allowed, please update alpha value and try again.'],alpha_max);
        uiwait(msgbox(message,'OK','modal'));
        handles.decon_success_flag = 0;
        return
    end
    save_folder = ['DeCon_Laguerre ' datestr(t,formatOut)];
elseif Decon_param.Decon_option == 3
    switch Decon_param.Exp_param.model
        case 1
            save_folder = ['DeCon_MonoExponential ' datestr(t,formatOut)];
            
        case 2
            save_folder = ['DeCon_BiExponential ' datestr(t,formatOut)];
            
        case 3
            save_folder = ['DeCon_TriExponential ' datestr(t,formatOut)];
    end
elseif Decon_param.Decon_option == 5
       save_folder = ['DeCon_AMD ' datestr(t,formatOut)];
else    
    warndlg('Deconvolution method not supported. Please double check your function input.','DeCon Method Warning.');
end

mkdir(folder_cell{1},save_folder); % make folder to save data
save_path = fullfile(folder_cell{1},save_folder); % set saving pass
%% initialize display cell array
raw_data_struct_cell = cell(size(name_cell));

INT1_map  = cell(size(name_cell));
INT1_R_map  = cell(size(name_cell));
LT1_map  = cell(size(name_cell));
INT2_map  = cell(size(name_cell));
INT2_R_map  = cell(size(name_cell));
LT2_map  = cell(size(name_cell));
INT3_map  = cell(size(name_cell));
INT3_R_map  = cell(size(name_cell));
LT3_map  = cell(size(name_cell));
INT4_map  = cell(size(name_cell));
INT4_R_map  = cell(size(name_cell));
LT4_map  = cell(size(name_cell));

INT1_map_interped  = cell(size(name_cell));
INT1_R_map_interped  = cell(size(name_cell));
LT1_map_interped  = cell(size(name_cell));
INT2_map_interped  = cell(size(name_cell));
INT2_R_map_interped  = cell(size(name_cell));
LT2_map_interped  = cell(size(name_cell));
INT3_map_interped  = cell(size(name_cell));
INT3_R_map_interped  = cell(size(name_cell));
LT3_map_interped  = cell(size(name_cell));
INT4_map_interped  = cell(size(name_cell));
INT4_R_map_interped  = cell(size(name_cell));
LT4_map_interped  = cell(size(name_cell));

SNR1_map  = cell(size(name_cell));
SNR2_map  = cell(size(name_cell));
SNR3_map  = cell(size(name_cell));
SNR4_map  = cell(size(name_cell));

SNR1_map_interped  = cell(size(name_cell));
SNR2_map_interped  = cell(size(name_cell));
SNR3_map_interped  = cell(size(name_cell));
SNR4_map_interped  = cell(size(name_cell));

CH1_Taus_1_map = cell(size(name_cell));
CH2_Taus_1_map = cell(size(name_cell));
CH3_Taus_1_map = cell(size(name_cell));
CH4_Taus_1_map = cell(size(name_cell));

CH1_Taus_2_map = cell(size(name_cell));
CH2_Taus_2_map = cell(size(name_cell));
CH3_Taus_2_map = cell(size(name_cell));
CH4_Taus_2_map = cell(size(name_cell));

CH1_Taus_1_map_interped = cell(size(name_cell));
CH2_Taus_1_map_interped = cell(size(name_cell));
CH3_Taus_1_map_interped = cell(size(name_cell));
CH4_Taus_1_map_interped = cell(size(name_cell));

CH1_Taus_2_map_interped = cell(size(name_cell));
CH2_Taus_2_map_interped = cell(size(name_cell));
CH3_Taus_2_map_interped = cell(size(name_cell));
CH4_Taus_2_map_interped = cell(size(name_cell));

CH1_Weight_1_map = cell(size(name_cell));
CH2_Weight_1_map = cell(size(name_cell));
CH3_Weight_1_map = cell(size(name_cell));
CH4_Weight_1_map = cell(size(name_cell));

CH1_Weight_2_map = cell(size(name_cell));
CH2_Weight_2_map = cell(size(name_cell));
CH3_Weight_2_map = cell(size(name_cell));
CH4_Weight_2_map = cell(size(name_cell));

CH1_Weight_1_map_interped = cell(size(name_cell));
CH2_Weight_1_map_interped = cell(size(name_cell));
CH3_Weight_1_map_interped = cell(size(name_cell));
CH4_Weight_1_map_interped = cell(size(name_cell));

CH1_Weight_2_map_interped = cell(size(name_cell));
CH2_Weight_2_map_interped = cell(size(name_cell));
CH3_Weight_2_map_interped = cell(size(name_cell));
CH4_Weight_2_map_interped = cell(size(name_cell));
% main GUI variable to save all the data structure
all_data = cell(size(name_cell,1),4);

% set progress bar to 0
set(handles.hPb,'Value',0);
global_shift = [];

%% added in function to align iRF based on global data rather than each line
if ~handles.irfAlignFlag
    temp_text = 'Global iRF alignment in process.';
    set(handles.edit_decon_process,'String',temp_text);
    set(handles.hPb,'Value',0); % update progress bar
    drawnow
    [handles, global_shift] = irfAlignGlobal(Decon_param, name_cell,folder_cell,Background,iIRF,handles);
    cd(save_path)
    %save all irf alignment result
    fileID = fopen('iRF_Alignment_all.txt','w');
    fprintf(fileID,'%6s %6s %6s %6s\r\n','CH1','CH2','CH3','CH4');
    fprintf(fileID,'%6.0f %6.0f %6.0f %6.0f\r\n',global_shift');
    fclose(fileID);
    global_shift = mode(global_shift,1);
    %save final irf align result
    fileID = fopen('iRF_Alignment.txt','w');
    fprintf(fileID,'%6s %6s %6s %6s\r\n','CH1','CH2','CH3','CH4');
    fprintf(fileID,'%6.0f %6.0f %6.0f %6.0f\r\n',global_shift);
    fclose(fileID);
%     set(handles.ch1Shift, 'enable', 'on');
%     set(handles.ch2Shift, 'enable', 'on');
%     set(handles.ch3Shift, 'enable', 'on');
%     set(handles.ch4Shift, 'enable', 'on');

else
    ch1shift = str2double(get(handles.ch1Shift, 'String'));
    ch2shift = str2double(get(handles.ch2Shift, 'String'));
    ch3shift = str2double(get(handles.ch3Shift, 'String'));
    ch4shift = str2double(get(handles.ch4Shift, 'String'));
    
    global_shift(1) = ch1shift;
    global_shift(2) = ch2shift;
    global_shift(3) = ch3shift;
    global_shift(4) = ch4shift;
    cd(save_path)
    fileID = fopen('iRF_Alignment.txt','w');
    fprintf(fileID,'%6s %6s %6s %6s\r\n','CH1','CH2','CH3','CH4');
    fprintf(fileID,'%6.0f %6.0f %6.0f %6.0f\r\n',global_shift');
    fclose(fileID);
end

%% loop for each line
for i=1:num_of_batches
    
    %% initilize all local variables
    % initilize common variables
    lifet_avg = cell(4,1);  % average lifetime
    lifet_avg_interped = cell(4,1); % interpreted average lifetime
    spec_int = cell(4,1);   % integrated intensity
    spec_int_R = cell(4,1);
    spec_int_interped = cell(4,1);   % integrated intensity
    spec_int_R_interped = cell(4,1);
    SNR = cell(4,1);        % SNR
    SNR_interped = cell(4,1);        % SNR
    test = cell(4,1);       % statistic structure
    Data = cell(4,1);       % cell array to save raw data structure
    Taus_1_cell = cell(4,1);
    Taus_2_cell = cell(4,1);
    Weights_1_cell = cell(4,1);
    Weights_2_cell = cell(4,1);
    Taus_1_cell_interped = cell(4,1);
    Taus_2_cell_interped = cell(4,1);
    Weights_1_cell_interped = cell(4,1);
    Weights_2_cell_interped = cell(4,1);
    % initilize method dependent variables
    if Decon_param.Decon_option == 1
        alpha_out = cell(4,1);  % alpha value for each channel
        Laguerre_coeffs = cell(4,1);    % Laguerre coefficinets
    elseif Decon_param.Decon_option == 3
        Taus_cell = cell(4,1);
        Weights_cell = cell(4,1);
        Taus_cell_interped = cell(4,1);
        Weights_cell_interped = cell(4,1);
    end
    %% preprocessing data
    
    rawDataClass = RawDataStrcutClass(name_cell(i),folder_cell(i),Background.data,iIRF.data); % creat raw data structure class
    n_pix = rawDataClass.rawdata_dimention(1); % get number of pixels per line
    ID = get(rawDataClass,'ID'); % ID of measurement
    timeStamp = get(rawDataClass,'timeStamp'); % get time stamp for interpoltion of the data
    spaceStamp = get(rawDataClass,'spaceStamp');% get space stmap for interpolation of data
    
    %     DAQtime = round(mean(diff(timeStamp))*1000); % calculate average DAQ time
    if all(~spaceStamp)
        Time = timeStamp(1)+(0:n_pix-1)*Decon_param.DAQ_time/1000; % acq time points for a frame
        Time = Time'; % reshape time
    else
        Time = spaceStamp;
    end
    peak_position = iIRF.peaks.*Decon_param.Channels;
    
    [pathstr,name,ext] = fileparts(name_cell{i}); % get file info
    
    if Decon_param.Decon_option == 1
        name = strcat(name_cell{i},'_Laguerre_DeCon.mat'); % creat saving name
    elseif Decon_param.Decon_option == 3
        switch Decon_param.Exp_param.model
            case 1
                %             savepath = fullfile(folder_cell{i},'DeCon_MonoExp');
                name = strcat(name_cell{i},'_MonoExp_DeCon.mat');
            case 2
                %             savepath = fullfile(folder_cell{i},'DeCon_BiExp');
                name = strcat(name_cell{i},'_BiExp_DeCon.mat');
            case 3
                %             savepath = fullfile(folder_cell{i},'DeCon_TriExp');
                name = strcat(name_cell{i},'_TriExp_DeCon.mat');
        end
    elseif Decon_param.Decon_option == 5
           name = strcat(name_cell{i},'AMD.mat'); % creat saving name 
    end
    % DC removal
    DC_range = getappdata(handles.figure_main,'DC_range');
    [~] = removeDC(rawDataClass,DC_range);
    % get gain information
    [gain_factor, gain_list] = gainCorrection(rawDataClass);
    
    [~,non_decon_idx,full_data_size]= detectSaturation(rawDataClass,Decon_param.amplitude_window(1),Decon_param.amplitude_window(2));
    
    
    decon_idx = rawDataClass.decon_idx; % get deconvolution index
    noise = get_noise(rawDataClass); %get noise form data
    
    %% set GUI progress bar
    temp_text = sprintf('Deconvolving Line %d with %d points',i,length(decon_idx));
    set(handles.edit_decon_process,'String',temp_text);
    set(handles.hPb,'Value',(i)/num_of_batches*100);
    drawnow
    
    %% initialize data matrix
    
    LT = zeros(full_data_size(1),1);
    INT = zeros(full_data_size(1),1);
    INT_R = zeros(full_data_size(1),1);
    BG_scale_factor = zeros(full_data_size(1),1);
    
    
    if Decon_param.Decon_option == 1
        CH_LCs = zeros(full_data_size(1),12);
        
        
    elseif Decon_param.Decon_option == 3
        Taus = zeros(full_data_size(1),exp_base_number);
        Weights = zeros(full_data_size(1),exp_base_number);
    end
    % check whether there is data to be deconvolve
    if ~isempty(decon_idx)
        
        % curser_low = 1000;
        % curser_high = 1200;
        
        
        if ~isempty(Background.refRange)
            [BG_removed, BG_scale_factor] = removeBG(rawDataClass,Background.refRange(1),Background.refRange(2));
            
        end
        
        if show_waveform_flag
            plot(BG_removed');
            ylim([-0.1 0.75])
            drawnow
        end
        
        
        datastruct = truncateData(rawDataClass,Decon_param.ch_width,peak_position); %cell array to store truncated data
        raw_data_struct_cell{i} = rawDataClass;
        
        if length(timeStamp) > full_data_size(1)
            timeStamp = timeStamp(1:full_data_size(1));
            Time = Time(1:full_data_size(1));
            gain_factor = gain_factor(1:full_data_size(1));
        end
        if length(timeStamp) < full_data_size(1)
            timeStamp = [0; timeStamp];
            Time = [0; Time];
            gain_factor = [1; gain_factor];
        end
        
        %check whether channel data is empty
        % DeCon loop
        for j = 1:4
            if ~isempty(datastruct{j})
                rawdatastruct = channeldata(datastruct{j}',iIRF.data{j},Decon_param.time_res,Decon_param.bw,decon_idx,noise,gain_list(decon_idx,:));
                SNR_temp = zeros(full_data_size(1),1);
                SNR_temp(decon_idx,:) = rawdatastruct.SNR;
                SNR{j} = SNR_temp;
                [~, index] = unique(timeStamp); 
                idx_temp = setdiff(1:size(timeStamp),index);
                timeStamp(idx_temp) = timeStamp(idx_temp)+eps;
                SNR_temp = interp1(timeStamp,SNR_temp,Time);
                SNR_interped{j} = SNR_temp;
                LT = zeros(full_data_size(1),1);
                %                 LT_interped = zeros(full_data_size(1),1);
                INT = zeros(full_data_size(1),1);
                INT_R = zeros(full_data_size(1),1);
                %                 INT_interped = zeros(full_data_size(1),1);
                if Decon_param.Decon_option == 1 % if user chose Laguerre
                    Laguerre_Struct=LaguerreModel(rawdatastruct,Decon_param.LG_param.alpha);
                    Laguerre_Struct.iIRF_align(mode(global_shift(:,j)));
%                     Laguerre_Struct.iIRF_align();
                    Laguerre_Struct.estimate_laguerre();
                    Data{j} = Laguerre_Struct;
                    all_data{i,j}=Laguerre_Struct;
                    
                    LT(decon_idx,:) = Laguerre_Struct.LTs';
                    LT_interped = interp1(timeStamp,LT,Time);
                    lifet_avg{j}=LT;
                    lifet_avg_interped{j}=LT_interped;
                    
                    % Weight1(decon_idx,:) = expsLaguerreModel_CH1.Weights';
                    INT(decon_idx,:) = Laguerre_Struct.INTs';
                    INT_R(decon_idx,:) = Laguerre_Struct.AUC';
                    INT = INT.*gain_factor;
                    INT_R = INT_R.*gain_factor;
                    INT_interped = interp1(timeStamp,INT,Time);
                    INT_R_interped = interp1(timeStamp,INT_R,Time);
                    spec_int{j}=INT;
                    spec_int_interped{j}=INT_interped;
                    spec_int_R{j} = INT_R;
                    spec_int_R_interped{j}=INT_R_interped;
                    CH_LCs(decon_idx,:) = Laguerre_Struct.LCs';
                    Laguerre_coeffs{j} = CH_LCs;
                    alpha_out{j} = get(Laguerre_Struct,'alpha');
                end
                if Decon_param.Decon_option == 5 % if user chose analog mean delay
%                     t1 = Decon_param.time_res*((1:size(rawdatastruct.data,1)))';
%                     t2 = Decon_param.time_res*((1:size(rawdatastruct.iIRF,1)))';
%                     AMD_data = rawdatastruct.data;
%                     AMD_iRF = circshift(rawdatastruct.iIRF,mode(global_shift(:,j)));
%                     
%                     INTs = sum(AMD_data,1);
%                     normalized_AMD_data = bsxfun(@rdivide,AMD_data,INTs)';
%                     delay_data = normalized_AMD_data*t1;
%                     delay_data(delay_data<0)=0;
%                     delay_iRF = AMD_iRF'*t2;
%                     LT_temp = delay_data-delay_iRF;
                    rawdatastruct.iIRF = circshift(rawdatastruct.iIRF,global_shift(j));
                    AMD_obj = AMDModel(rawdatastruct);
                    AMD_obj.calculateLT;
                    Data{j} = AMD_obj;
                    LT(decon_idx,:) = AMD_obj.LTs;
                    LT_interped = interp1(timeStamp,LT,Time);
                    lifet_avg{j}=LT;
                    lifet_avg_interped{j}=LT_interped;
                    
                    % Weight1(decon_idx,:) = expsLaguerreModel_CH1.Weights';
                    INT(decon_idx,:) = AMD_obj.INTs;
                    INT = INT.*gain_factor;
                    INT_interped = interp1(timeStamp,INT,Time);
                    spec_int{j}=INT;
                    spec_int_interped{j}=INT_interped;
                                       
                end
                if Decon_param.Decon_option == 3 % if user chose Exponential
                    Taus = zeros(full_data_size(1),1);
                    %                     Taus_interped = zeros(full_data_size(1),1);
                    Weights = zeros(full_data_size(1),1);
                    %                     Weights_interped = zeros(full_data_size(1),1);
                    ExpStruct=ExpModel(rawdatastruct,Decon_param.Exp_param.model);
                    ExpStruct.iIRF_align();
                    ExpStruct.fit_exp();
                    Data{j} = ExpStruct;
                    all_data{i,j}=ExpStruct;
                    
                    LT(decon_idx,:) = ExpStruct.LTs';
                    LT_interped = interp1(timeStamp,LT,Time);
                    lifet_avg{j}=LT;
                    lifet_avg_interped{j}=LT_interped;
                    % Weight1(decon_idx,:) = expsLaguerreModel_CH1.Weights';
                    INT(decon_idx,:) = ExpStruct.INTs';
                    INT = INT.*gain_factor;
                    INT_interped = interp1(timeStamp,INT,Time);
                    spec_int{j}=INT;
                    spec_int_interped{j}=INT_interped;
                    
                    Taus(decon_idx,1:exp_base_number) = ExpStruct.Taus';
                    Taus_interped = interp1(timeStamp,Taus,Time);
                    Taus_1_cell{j} = Taus(:,1);
                    Taus_1_cell_interped{j} = Taus_interped(:,1);
                    try
                        Taus_2_cell{j} = Taus(:,2);
                        Taus_2_cell_interped{j} = Taus_interped(:,2);
                    catch
                    end
                    Weights(decon_idx,1:exp_base_number) = ExpStruct.Weights';
                    Weights_interped = interp1(timeStamp,Weights,Time);
                    Weights_1_cell{j} = Weights(:,1);
                    Weights_1_cell_interped{j} = Weights_interped(:,1);
                    try
                        Weights_2_cell{j} = Weights(:,2);
                        Weights_2_cell_interped{j} = Weights_interped(:,2);
                    catch
                    end
                end
                % clear variables for each line
                clear SNR_temp LT LT_interped INT INT_interped CH_LCs LaguerreStruct Taus Taus_interped Weights Weights_interped
            end
            
        end
        
        
        
    else % if nothing to deconvolve, set everything to 0
        
        for k=1:4
            spec_int{k} = zeros(full_data_size(1),1);
            spec_int_R{k} = zeros(full_data_size(1),1);
            lifet_avg{k} = zeros(full_data_size(1),1);
            SNR{k} = zeros(full_data_size(1),1);
            spec_int_interped{k} = zeros(full_data_size(1),1);
            spec_int_R_interped = zeros(full_data_size(1),1);
            lifet_avg_interped{k} = zeros(full_data_size(1),1);
            SNR_interped{k} = zeros(full_data_size(1),1);
        end
    end
    %% save data to cell arrays
    INT1_map{i}  = spec_int{1}';
    INT1_R_map{i}  = spec_int_R{1}';
    LT1_map{i}  = lifet_avg{1}';
    SNR1_map{i}  = SNR{1}';
    INT2_map{i} = spec_int{2}';
    INT2_R_map{i} = spec_int_R{2}';
    LT2_map{i}  = lifet_avg{2}';
    SNR2_map{i}  = SNR{2}';
    INT3_map{i}  = spec_int{3}';
    INT3_R_map{i} = spec_int_R{3}';
    LT3_map{i}  = lifet_avg{3}';
    SNR3_map{i}  = SNR{3}';
    INT4_map{i} = spec_int{4}';
    INT4_R_map{i} = spec_int_R{4}';
    LT4_map{i}  = lifet_avg{4}';
    SNR4_map{i}  = SNR{4}';
    
    INT1_map_interped{i}  = spec_int_interped{1}';
    INT1_R_map_interped{i} = spec_int_R_interped{1}';
    LT1_map_interped{i}  = lifet_avg_interped{1}';
    SNR1_map_interped{i}  = SNR_interped{1}';
    INT2_map_interped{i} = spec_int_interped{2}';
    INT2_R_map_interped{i} = spec_int_R_interped{2}';
    LT2_map_interped{i}  = lifet_avg_interped{2}';
    SNR2_map_interped{i}  = SNR_interped{2}';
    INT3_map_interped{i}  = spec_int_interped{3}';
    INT3_R_map_interped{i} = spec_int_R_interped{3}';
    LT3_map_interped{i}  = lifet_avg_interped{3}';
    SNR3_map_interped{i}  = SNR_interped{3}';
    INT4_map_interped{i} = spec_int_interped{4}';
    INT4_R_map_interped{i} = spec_int_R_interped{4}';
    LT4_map_interped{i}  = lifet_avg_interped{4}';
    SNR4_map_interped{i}  = SNR_interped{4}';
    
    if Decon_param.Decon_option == 3
        
        %         Taus_1_cell = zeros(full_data_size(1),1);
        %         Weights_1_cell = zeros(full_data_size(1),1);
        %         Taus_1_cell_interped  = zeros(full_data_size(1),1);
        %         Weights_1_cell_interped = zeros(full_data_size(1),1);
        
        CH1_Taus_1_map{i} = Taus_1_cell{1}';
        CH2_Taus_1_map{i} = Taus_1_cell{2}';
        CH3_Taus_1_map{i} = Taus_1_cell{3}';
        CH4_Taus_1_map{i} = Taus_1_cell{4}';
        
        CH1_Weight_1_map{i} = Weights_1_cell{1}';
        CH2_Weight_1_map{i} = Weights_1_cell{2}';
        CH3_Weight_1_map{i} = Weights_1_cell{3}';
        CH4_Weight_1_map{i} = Weights_1_cell{4}';
        
        CH1_Taus_1_map_interped{i} = Taus_1_cell_interped{1}';
        CH2_Taus_1_map_interped{i} = Taus_1_cell_interped{2}';
        CH3_Taus_1_map_interped{i} = Taus_1_cell_interped{3}';
        CH4_Taus_1_map_interped{i} = Taus_1_cell_interped{4}';
        
        CH1_Weight_1_map_interped{i} = Weights_1_cell_interped{1}';
        CH2_Weight_1_map_interped{i} = Weights_1_cell_interped{2}';
        CH3_Weight_1_map_interped{i} = Weights_1_cell_interped{3}';
        CH4_Weight_1_map_interped{i} = Weights_1_cell_interped{4}';
        
        
        try
            CH1_Taus_2_map{i} = Taus_2_cell{1}';
            CH2_Taus_2_map{i} = Taus_2_cell{2}';
            CH3_Taus_2_map{i} = Taus_2_cell{3}';
            CH4_Taus_2_map{i} = Taus_2_cell{4}';
            
            CH1_Weight_2_map{i} = Weights_2_cell{1}';
            CH2_Weight_2_map{i} = Weights_2_cell{2}';
            CH3_Weight_2_map{i} = Weights_2_cell{3}';
            CH4_Weight_2_map{i} = Weights_2_cell{4}';
            
            CH1_Taus_2_map_interped{i} = Taus_2_cell_interped{1}';
            CH2_Taus_2_map_interped{i} = Taus_2_cell_interped{2}';
            CH3_Taus_2_map_interped{i} = Taus_2_cell_interped{3}';
            CH4_Taus_2_map_interped{i} = Taus_2_cell_interped{4}';
            
            CH1_Weight_2_map_interped{i} = Weights_2_cell_interped{1}';
            CH2_Weight_2_map_interped{i} = Weights_2_cell_interped{2}';
            CH3_Weight_2_map_interped{i} = Weights_2_cell_interped{3}';
            CH4_Weight_2_map_interped{i} = Weights_2_cell_interped{4}';
        catch
        end
    end
    
    
    
    %% save to .mat file
    cd(save_path)
    BG = Background.data;
    if Decon_param.Decon_option == 1
        save(name,'lifet_avg','lifet_avg_interped','spec_int','spec_int_interped',...
            'alpha_out','Laguerre_coeffs',...
            'timeStamp','SNR','SNR_interped','test','decon_idx','full_data_size','BG',...
            'BG_scale_factor','gain_factor','gain_list','Data','spaceStamp','ID','rawDataClass');
    end
    if Decon_param.Decon_option == 3
        save(name,'lifet_avg','lifet_avg_interped','spec_int','spec_int_interped',...
            'Weights_1_cell','Weights_1_cell_interped','Weights_2_cell','Weights_2_cell_interped',...
            'Taus_1_cell','Taus_1_cell_interped','Taus_2_cell','Taus_2_cell_interped',...
            'timeStamp','SNR','SNR_interped','test','decon_idx','full_data_size','BG',...
            'BG_scale_factor','gain_factor','gain_list','Data','spaceStamp','ID','rawDataClass');
    end
    if Decon_param.Decon_option == 5
        save(name,'lifet_avg','lifet_avg_interped','spec_int','spec_int_interped',...
            'timeStamp','SNR','SNR_interped','test','decon_idx','full_data_size','BG',...
            'BG_scale_factor','gain_factor','gain_list','Data','spaceStamp','ID','rawDataClass');
    end
    
end % end of the loop of line

%% pass data to main GUI

INT_map{1} = pad_zero_and_stitch(INT1_map);
INT_map{2} = pad_zero_and_stitch(INT2_map);
INT_map{3} = pad_zero_and_stitch(INT3_map);
INT_map{4} = pad_zero_and_stitch(INT4_map);
INT_R_map{1} = pad_zero_and_stitch(INT1_R_map);
INT_R_map{2} = pad_zero_and_stitch(INT2_R_map);
INT_R_map{3} = pad_zero_and_stitch(INT3_R_map);
INT_R_map{4} = pad_zero_and_stitch(INT4_R_map);
INT_map_interped{1} = pad_zero_and_stitch(INT1_map_interped);
INT_map_interped{2} = pad_zero_and_stitch(INT2_map_interped);
INT_map_interped{3} = pad_zero_and_stitch(INT3_map_interped);
INT_map_interped{4} = pad_zero_and_stitch(INT4_map_interped);
INT_R_map_interped{1} = pad_zero_and_stitch(INT1_R_map_interped);
INT_R_map_interped{2} = pad_zero_and_stitch(INT2_R_map_interped);
INT_R_map_interped{3} = pad_zero_and_stitch(INT3_R_map_interped);
INT_R_map_interped{4} = pad_zero_and_stitch(INT4_R_map_interped);
setappdata(handles.figure_main,'INT_map',INT_map);
setappdata(handles.figure_main,'INT_R_map',INT_R_map);
setappdata(handles.figure_main,'INT_map_interped',INT_map_interped);
setappdata(handles.figure_main,'INT_R_map_interped',INT_R_map_interped);

LT_map{1} = pad_zero_and_stitch(LT1_map);
LT_map{2} = pad_zero_and_stitch(LT2_map);
LT_map{3} = pad_zero_and_stitch(LT3_map);
LT_map{4} = pad_zero_and_stitch(LT4_map);
LT_map_interped{1} = pad_zero_and_stitch(LT1_map_interped);
LT_map_interped{2} = pad_zero_and_stitch(LT2_map_interped);
LT_map_interped{3} = pad_zero_and_stitch(LT3_map_interped);
LT_map_interped{4} = pad_zero_and_stitch(LT4_map_interped);
setappdata(handles.figure_main,'LT_map',LT_map);
setappdata(handles.figure_main,'LT_map_interped',LT_map_interped);

SNR_map{1} = pad_zero_and_stitch(SNR1_map);
SNR_map{2} = pad_zero_and_stitch(SNR2_map);
SNR_map{3} = pad_zero_and_stitch(SNR3_map);
SNR_map{4} = pad_zero_and_stitch(SNR4_map);
SNR_map_interped{1} = pad_zero_and_stitch(SNR1_map_interped);
SNR_map_interped{2} = pad_zero_and_stitch(SNR2_map_interped);
SNR_map_interped{3} = pad_zero_and_stitch(SNR3_map_interped);
SNR_map_interped{4} = pad_zero_and_stitch(SNR4_map_interped);

setappdata(handles.figure_main,'SNR_map',SNR_map);
setappdata(handles.figure_main,'SNR_map_interped',SNR_map_interped);

if Decon_param.Decon_option == 3
    Taus_1_map{1} = pad_zero_and_stitch(CH1_Taus_1_map);
    Taus_1_map{2} = pad_zero_and_stitch(CH2_Taus_1_map);
    Taus_1_map{3} = pad_zero_and_stitch(CH3_Taus_1_map);
    Taus_1_map{4} = pad_zero_and_stitch(CH4_Taus_1_map);
    setappdata(handles.figure1,'Taus_1_map',Taus_1_map);
    
    
    Taus_2_map{1} = pad_zero_and_stitch(CH1_Taus_2_map);
    Taus_2_map{2} = pad_zero_and_stitch(CH2_Taus_2_map);
    Taus_2_map{3} = pad_zero_and_stitch(CH3_Taus_2_map);
    Taus_2_map{4} = pad_zero_and_stitch(CH4_Taus_2_map);
    setappdata(handles.figure1,'Taus_2_map',Taus_2_map);
    
    
    Taus_3_map{1} = pad_zero_and_stitch(CH1_Taus_3_map);
    Taus_3_map{2} = pad_zero_and_stitch(CH2_Taus_3_map);
    Taus_3_map{3} = pad_zero_and_stitch(CH3_Taus_3_map);
    Taus_3_map{4} = pad_zero_and_stitch(CH4_Taus_3_map);
    setappdata(handles.figure1,'Taus_3_map',Taus_3_map);
    
    Taus_1_map_interped{1} = pad_zero_and_stitch(CH1_Taus_1_map_interped);
    Taus_1_map_interped{2} = pad_zero_and_stitch(CH2_Taus_1_map_interped);
    Taus_1_map_interped{3} = pad_zero_and_stitch(CH3_Taus_1_map_interped);
    Taus_1_map_interped{4} = pad_zero_and_stitch(CH4_Taus_1_map_interped);
    setappdata(handles.figure1,'Taus_1_map_interped',Taus_1_map_interped);
    
    Taus_2_map_interped{1} = pad_zero_and_stitch(CH1_Taus_2_map_interped);
    Taus_2_map_interped{2} = pad_zero_and_stitch(CH2_Taus_2_map_interped);
    Taus_2_map_interped{3} = pad_zero_and_stitch(CH3_Taus_2_map_interped);
    Taus_2_map_interped{4} = pad_zero_and_stitch(CH4_Taus_2_map_interped);
    setappdata(handles.figure1,'Taus_2_map_interped',Taus_2_map_interped);
    
    Taus_3_map_interped{1} = pad_zero_and_stitch(CH1_Taus_3_map_interped);
    Taus_3_map_interped{2} = pad_zero_and_stitch(CH2_Taus_3_map_interped);
    Taus_3_map_interped{3} = pad_zero_and_stitch(CH3_Taus_3_map_interped);
    Taus_3_map_interped{4} = pad_zero_and_stitch(CH4_Taus_3_map_interped);
    setappdata(handles.figure1,'Taus_3_map_interped',Taus_3_map_interped);
    
    %save weights map
    Weight_1_map{1} = pad_zero_and_stitch(CH1_Weight_1_map);
    Weight_1_map{2} = pad_zero_and_stitch(CH2_Weight_1_map);
    Weight_1_map{3} = pad_zero_and_stitch(CH3_Weight_1_map);
    Weight_1_map{4} = pad_zero_and_stitch(CH4_Weight_1_map);
    setappdata(handles.figure1,'Weight_1_map',Weight_1_map);
    
    Weight_2_map{1} = pad_zero_and_stitch(CH1_Weight_2_map);
    Weight_2_map{2} = pad_zero_and_stitch(CH2_Weight_2_map);
    Weight_2_map{3} = pad_zero_and_stitch(CH3_Weight_2_map);
    Weight_2_map{4} = pad_zero_and_stitch(CH4_Weight_2_map);
    setappdata(handles.figure1,'Weight_2_map',Weight_2_map);
    
    Weight_3_map{1} = pad_zero_and_stitch(CH1_Weight_3_map);
    Weight_3_map{2} = pad_zero_and_stitch(CH2_Weight_3_map);
    Weight_3_map{3} = pad_zero_and_stitch(CH3_Weight_3_map);
    Weight_3_map{4} = pad_zero_and_stitch(CH4_Weight_3_map);
    setappdata(handles.figure1,'Weight_3_map',Weight_3_map);
    
    Weight_1_map_interped{1} = pad_zero_and_stitch(CH1_Weight_1_map_interped);
    Weight_1_map_interped{2} = pad_zero_and_stitch(CH2_Weight_1_map_interped);
    Weight_1_map_interped{3} = pad_zero_and_stitch(CH3_Weight_1_map_interped);
    Weight_1_map_interped{4} = pad_zero_and_stitch(CH4_Weight_1_map_interped);
    setappdata(handles.figure1,'Weight_1_map_interped',Weight_1_map_interped);
    
    Weight_2_map_interped{1} = pad_zero_and_stitch(CH1_Weight_2_map_interped);
    Weight_2_map_interped{2} = pad_zero_and_stitch(CH2_Weight_2_map_interped);
    Weight_2_map_interped{3} = pad_zero_and_stitch(CH3_Weight_2_map_interped);
    Weight_2_map_interped{4} = pad_zero_and_stitch(CH4_Weight_2_map_interped);
    setappdata(handles.figure1,'Weight_2_map_interped',Weight_2_map_interped);
    
    Weight_3_map_interped{1} = pad_zero_and_stitch(CH1_Weight_3_map_interped);
    Weight_3_map_interped{2} = pad_zero_and_stitch(CH2_Weight_3_map_interped);
    Weight_3_map_interped{3} = pad_zero_and_stitch(CH3_Weight_3_map_interped);
    Weight_3_map_interped{4} = pad_zero_and_stitch(CH4_Weight_3_map_interped);
    
    setappdata(handles.figure1,'Weight_3_map_interped',Weight_3_map_interped);
    
    
    setappdata(handles.figure_main,'Taus_1_map',Taus_1_map);
    setappdata(handles.figure_main,'Taus_2_map',Taus_2_map);
    setappdata(handles.figure_main,'Taus_3_map',Taus_2_map);
    
    setappdata(handles.figure_main,'Taus_1_map_interped',Taus_1_map_interped);
    setappdata(handles.figure_main,'Taus_2_map_interped',Taus_2_map_interped);
    setappdata(handles.figure_main,'Taus_3_map_interped',Taus_1_map_interped);
    setappdata(handles.figure_main,'Weight_1_map',Weight_1_map);
    setappdata(handles.figure_main,'Weight_2_map',Weight_2_map);
    setappdata(handles.figure_main,'Weight_3_map',Weight_2_map);
    
    setappdata(handles.figure_main,'Weight_1_map_interped',Weight_1_map_interped);
    setappdata(handles.figure_main,'Weight_2_map_interped',Weight_2_map_interped);
    setappdata(handles.figure_main,'Weight_3_map_interped',Weight_1_map_interped);
    
    
    
    
    save([name_cell{1}(1:end-3) '_reconstructed'],'INT_map','LT_map','SNR_map',...
        'INT_map_interped','LT_map_interped','SNR_map_interped','Taus_1_map',...
        'Taus_2_map','Taus_1_map_interped','Taus_2_map_interped',...
        'Taus_3_map','Taus_3_map_interped','Weight_1_map','Weight_2_map',...
        'Weight_1_map_interped','Weight_2_map_interped','Weight_3_map',...
        'Weight_3_map_interped');
else
    save([name_cell{1}(1:end-3) '_reconstructed'],'INT_map','LT_map','SNR_map',...
        'INT_map_interped','LT_map_interped','SNR_map_interped','INT_R_map','INT_R_map_interped');
end


setappdata(handles.figure_main,'all_data',all_data);
setappdata(handles.figure_main,'raw_data_cell',raw_data_struct_cell);
handles.irfAlignFlag = 1;
handles.decon_success_flag = 1;
cd(handles.rootFolder)

