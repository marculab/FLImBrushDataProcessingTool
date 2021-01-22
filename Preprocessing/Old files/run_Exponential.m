function handles = run_Exponential(name_cell, folder_cell,BG,irf,handles,show_waveform_flag)
% get necessary parameter
low = getappdata(handles.figure1,'low');
high = getappdata(handles.figure1,'high');
dt = getappdata(handles.figure1,'dt');
BW = getappdata(handles.figure1,'BW');
truncation_length = getappdata(handles.figure1,'truncation_length');
allItems = get(handles.popupmenu7,'string');
selectedIndex = get(handles.popupmenu7,'Value');
exp_base_number = str2double(allItems{selectedIndex});
% make saving directery
num_of_batches = length(name_cell);
% get current date and time
t=datetime('now');
formatOut = 'mmm-dd-yyyy HH.MM PM';

switch exp_base_number
    case 1
        save_folder = ['DeCon_MomoExponential ' datestr(t,formatOut)];
        
    case 2
        save_folder = ['DeCon_BiExponential ' datestr(t,formatOut)];
        
    case 3
        save_folder = ['DeCon_TriExponential ' datestr(t,formatOut)];
        
end
% make new directory
mkdir(folder_cell{1},save_folder);
savepath = fullfile(folder_cell{1},save_folder);

%initialize diesplay cell array
INT1_map  = cell(size(name_cell));
LT1_map  = cell(size(name_cell));
INT2_map  = cell(size(name_cell));
LT2_map  = cell(size(name_cell));
INT3_map  = cell(size(name_cell));
LT3_map  = cell(size(name_cell));
INT4_map  = cell(size(name_cell));
LT4_map  = cell(size(name_cell));

SNR1_map  = cell(size(name_cell));
SNR2_map  = cell(size(name_cell));
SNR3_map  = cell(size(name_cell));
SNR4_map  = cell(size(name_cell));


handles.all_data = cell(size(name_cell,1),4);
% specify ploting axis
axes(handles.ax_decay);
set(handles.hPb,'Value',0);
%% loop for each line
for i=1:num_of_batches
    
    temp_text = sprintf('Deconvolving Line %d',i);
    set(handles.Decon_process_text,'String',temp_text);
    drawnow
    
    
    lifet_avg = cell(4,1);
    spec_int = cell(4,1);
    %     alpha_out = cell(4,1);
    SNR = cell(4,1);
    test = cell(4,1);
    Taus = cell(4,1);
    Weights = cell(4,1);
    Data = cell(4,1);
    data = RawDataStrcutClass(name_cell(i),folder_cell(i),BG,irf);
    
    BG_remove_postion = getappdata(handles.figure1,'BG_remove_postion');
    peak_position = getappdata(handles.figure1,'peak_position');
    
    
    timeStamp = get(data,'timeStamp');
    [pathstr,name,ext] = fileparts(name_cell{i});
    %     savepath = fullfile(folder_cell{i},'DeCon_Laguerre',strcat(name_cell{i},'_Laguerre_DeCon.mat'));
    % get gain information
    switch exp_base_number
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
    
    [gain_factor, gain_list] = gainCorrection(data);
    
    [~,non_decon_idx,full_data_size]= detectSaturation(data,low,high);
    
    decon_idx = data.decon_idx;
    % reduced_data = out;
    % reduced_data(idx,:)=[];
    % currentdata  = get(data,'rawdata');
    % dif = currentdata-out;
    
    %get noise from the tail of the waveform
    noise = get_noise(data);
    
    %initialize final data matrix
    LT1 = zeros(full_data_size(1),1);
    LT2 = zeros(full_data_size(1),1);
    LT3 = zeros(full_data_size(1),1);
    LT4 = zeros(full_data_size(1),1);
    % Weight1 = zeros(full_data_size(1),2);
    % Weight2 = zeros(full_data_size(1),2);
    INT1 = zeros(full_data_size(1),1);
    INT2 = zeros(full_data_size(1),1);
    INT3 = zeros(full_data_size(1),1);
    INT4 = zeros(full_data_size(1),1);
    
    
    BG_scale_factor = zeros(full_data_size(1),1);
    % Tau1 = zeros(full_data_size(1),2);
    % Tau2 = zeros(full_data_size(1),2);
    CH1_Taus = zeros(full_data_size(1),exp_base_number);
    CH2_Taus = zeros(full_data_size(1),exp_base_number);
    CH3_Taus = zeros(full_data_size(1),exp_base_number);
    CH4_Taus = zeros(full_data_size(1),exp_base_number);
    
    CH1_Weights = zeros(full_data_size(1),exp_base_number);
    CH2_Weights = zeros(full_data_size(1),exp_base_number);
    CH3_Weights = zeros(full_data_size(1),exp_base_number);
    CH4_Weights = zeros(full_data_size(1),exp_base_number);
    
    
    % check whether there is data to be deconvolve
    if ~isempty(decon_idx)
        
        % curser_low = 1000;
        % curser_high = 1200;
        
        [BG_removed, BG_scale_factor] = removeBG(data,BG_remove_postion(1),BG_remove_postion(2));
        
        if show_waveform_flag
            plot(BG_removed');
            ylim([-0.1 0.75])
            drawnow
        end
        % handles.truncation_length=680;
        % peak_position =[515 1281 0 0];
        
        datastruct = truncateData(data,truncation_length,peak_position);
        
        % dt = 0.08;
        
        % data = datasruct{1}';
        
        % mask = max(datasruct{1},[],2);
        % mask=mask>0.05;
        % data = data(:,mask);
        % [A, T, avglife, intensity, fitt, raw] = multiexp_fit(data,dt,laser,2);
        % addpath(genpath('C:\Users\XZ\OneDrive\UC Davis\Projects\Programming\Matlab Code\trfs_triplex_gui'));
        
        %check whether channel data is empty
        %% CH1 deconvolution
        if ~isempty(datastruct{1})
            rawdatatruct1 = channeldata(datastruct{1}',irf{1},dt,BW,decon_idx,noise,gain_list);
            SNR_temp = zeros(full_data_size(1),1);
            SNR_temp(decon_idx,:) = rawdatatruct1.SNR;
            SNR{1} = SNR_temp;
            clear SNR_temp
            ExpStruct_CH1=ExpModel(rawdatatruct1,exp_base_number);
            ExpStruct_CH1.iIRF_align();
            ExpStruct_CH1.fit_exp();
            
            handles.all_data{i,1}=ExpStruct_CH1;
            
            LT1(decon_idx,:) = ExpStruct_CH1.LTs';
            lifet_avg{1}=LT1;
            % Weight1(decon_idx,:) = expsLaguerreModel_CH1.Weights';
            INT1(decon_idx,:) = ExpStruct_CH1.INTs';
            INT1 = INT1.*gain_factor;
            spec_int{1}=INT1;
            CH1_Taus(decon_idx,1:exp_base_number) = ExpStruct_CH1.Taus';
            Taus{1} = CH1_Taus;
            CH1_Weights(decon_idx,1:exp_base_number) = ExpStruct_CH1.Weights';
            Weights{1} = CH1_Weights;
        end
        % h=figure;
        % plot(rawdatatruct1.data);
        % hold on
        %
        % plot(get(LaguerreModel_CH1,'fit'));
        
        %% CH2 deconvolution
        if ~isempty(datastruct{2})
            rawdatatruct2 = channeldata(datastruct{2}',irf{2},dt,BW,decon_idx,noise,gain_list);
            SNR_temp = zeros(full_data_size(1),1);
            SNR_temp(decon_idx,:) = rawdatatruct2.SNR;
            SNR{2} = SNR_temp;
            clear SNR_temp
            ExpStruct_CH2=ExpModel(rawdatatruct2,exp_base_number);
            ExpStruct_CH2.iIRF_align();
            ExpStruct_CH2.fit_exp();
            
            handles.all_data{i,2}=ExpStruct_CH2;
            
            LT2(decon_idx,:) = ExpStruct_CH2.LTs';
            lifet_avg{2}=LT2;

            % Weight1(decon_idx,:) = expsLaguerreModel_CH1.Weights';
            INT2(decon_idx,:) = ExpStruct_CH2.INTs';
            INT2 = INT2.*gain_factor;
            spec_int{2}=INT2;
            CH2_Taus(decon_idx,1:exp_base_number) = ExpStruct_CH2.Taus';
            Taus{2} = CH2_Taus;
            CH2_Weights(decon_idx,1:exp_base_number) = ExpStruct_CH2.Weights';
            Weights{2} = CH2_Weights;
        end
        
        %% CH3 deconvolution
        
        if ~isempty(datastruct{3})
            rawdatatruct3 = channeldata(datastruct{3}',irf{3},dt,BW,decon_idx,noise,gain_list);
            SNR_temp = zeros(full_data_size(1),1);
            SNR_temp(decon_idx,:) = rawdatatruct3.SNR;
            SNR{3} = SNR_temp;
            clear SNR_temp
            ExpStruct_CH3=ExpModel(rawdatatruct3,exp_base_number);
            ExpStruct_CH3.iIRF_align();
            ExpStruct_CH3.fit_exp();
            
            handles.all_data{i,3}=ExpStruct_CH3;
            
            LT3(decon_idx,:) = ExpStruct_CH3.LTs';
            lifet_avg{3}=LT3;
            % Weight1(decon_idx,:) = expsLaguerreModel_CH1.Weights';
            INT3(decon_idx,:) = ExpStruct_CH3.INTs';
            INT3 = INT3.*gain_factor;
            spec_int{3}=INT3;
            CH3_Taus(decon_idx,1:exp_base_number) = ExpStruct_CH3.Taus';
            Taus{3} = CH3_Taus;
            CH3_Weights(decon_idx,1:exp_base_number) = ExpStruct_CH3.Weights';
            Weights{3} = CH3_Weights;
        end
        
        
        %% CH4 deconvolution
        if ~isempty(datastruct{4})
            rawdatatruct4 = channeldata(datastruct{4}',irf{4},dt,BW,decon_idx,noise,gain_list);
            SNR_temp = zeros(full_data_size(1),1);
            SNR_temp(decon_idx,:) = rawdatatruct4.SNR;
            SNR{4} = SNR_temp;
            clear SNR_temp
            ExpStruct_CH4=ExpModel(rawdatatruct4,exp_base_number);
            ExpStruct_CH4.iIRF_align();
            ExpStruct_CH4.fit_exp();
            
            handles.all_data{i,4}=ExpStruct_CH4;
            
            LT4(decon_idx,:) = ExpStruct_CH4.LTs';
            lifet_avg{4}=LT4;
            % Weight1(decon_idx,:) = expsLaguerreModel_CH1.Weights';
            INT4(decon_idx,:) = ExpStruct_CH4.INTs';
            INT4 = INT4.*gain_factor;
            spec_int{4}=INT4;
            CH4_Taus(decon_idx,1:exp_base_number) = ExpStruct_CH4.Taus';
            Taus{4} = CH4_Taus;
            CH4_Weights(decon_idx,1:exp_base_number) = ExpStruct_CH4.Weights';
            Weights{4} = CH4_Weights;
        end
        
        
        % Tau1(decon_idx,:) = expsLaguerreModel_CH1.Taus';
        % Tau2(decon_idx,:) = expstruct2.Taus';
        
    end
    
    
    %cell array to save data class
    if exist('LaguerreModel_CH1','var')
    Data{1} = LaguerreModel_CH1;
    end
    if exist('LaguerreModel_CH2','var')
    Data{2} = LaguerreModel_CH2;
    end
    if exist('LaguerreModel_CH3','var')
    Data{3} = LaguerreModel_CH3;
    end
    if exist('LaguerreModel_CH4','var')
    Data{4} = LaguerreModel_CH4;
    end
    
    % pass data to main GUI
    INT1_map{i}  = INT1';
    LT1_map{i}  = LT1';
    SNR1_map{i}  = SNR{1}';
    INT2_map{i} = INT2';
    LT2_map{i}  = LT2';
    SNR2_map{i}  = SNR{2}';
    INT3_map{i}  = INT3';
    LT3_map{i}  = LT3';
    SNR3_map{i}  = SNR{3}';
    INT4_map{i} = INT4';
    LT4_map{i}  = LT4';
    SNR4_map{i}  = SNR{4}';
    
    handles.INT_map{1} = cell2mat(INT1_map);
    handles.INT_map{2} = cell2mat(INT2_map);
    handles.INT_map{3} = cell2mat(INT3_map);
    handles.INT_map{4} = cell2mat(INT4_map);
    
    handles.LT_map{1} = cell2mat(LT1_map);
    handles.LT_map{2} = cell2mat(LT2_map);
    handles.LT_map{3} = cell2mat(LT3_map);
    handles.LT_map{4} = cell2mat(LT4_map);
    
    handles.SNR_map{1} = cell2mat(SNR1_map);
    handles.SNR_map{2} = cell2mat(SNR2_map);
    handles.SNR_map{3} = cell2mat(SNR3_map);
    handles.SNR_map{4} = cell2mat(SNR4_map);
    
    
    % save to .mat file
    cd(savepath)
    save(name,'lifet_avg','spec_int','Weights','Taus',...
        'timeStamp','SNR','test','decon_idx','full_data_size','BG',...
        'BG_scale_factor','gain_factor','gain_list','Data');
    
    set(handles.hPb,'Value',i/num_of_batches*100);
    
end
    cd(handles.rootFolder)