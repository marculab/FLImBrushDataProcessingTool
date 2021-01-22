function handles = run_Laguerre(name_cell, folder_cell,BG,irf,handles,show_waveform_flag)
% get necessary parameter
low = getappdata(handles.figure1,'low');
high = getappdata(handles.figure1,'high');
dt = getappdata(handles.figure1,'dt');
BW = getappdata(handles.figure1,'BW');
truncation_length = getappdata(handles.figure1,'truncation_length');
% make saving directery
num_of_batches = length(name_cell);
t=datetime('now');
formatOut = 'mmm-dd-yyyy HH.MM PM';
save_folder = ['DeCon_Laguerre ' datestr(t,formatOut)];
mkdir(folder_cell{1},save_folder);

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
% main GUI variable to save all the data structure
handles.all_data = cell(size(name_cell,1),4);
% specify ploting axis
axes(handles.ax_decay);
% set progress bar to 0
set(handles.hPb,'Value',0);
%% loop for each line
for i=1:num_of_batches
    % set GUI text
    
    % initilize all local variables
    lifet_avg = cell(4,1);  % average lifetime
    spec_int = cell(4,1);   % integrated intensity
    alpha_out = cell(4,1);  % alpha value for each channel
    SNR = cell(4,1);        % SNR
    test = cell(4,1);       % statistic structure
    Laguerre_coeffs = cell(4,1);    % Laguerre coefficinets
    Data = cell(4,1);       % cell array to save raw data structure
    
    data = RawDataStrcutClass(name_cell(i),folder_cell(i),BG,irf);
    n_pix = data.rawdata_dimention(1);
    timeStamp = get(data,'timeStamp');
    DAQtime = round(mean(diff(timeStamp))*1000);
    Time = (0:n_pix-1)*DAQtime/1000; % acq time points for a frame
    Time = Time';
    BG_remove_postion = getappdata(handles.figure1,'BG_remove_postion');
    peak_position = getappdata(handles.figure1,'peak_position');

    
    [pathstr,name,ext] = fileparts(name_cell{i});
    %     savepath = fullfile(folder_cell{i},'DeCon_Laguerre',strcat(name_cell{i},'_Laguerre_DeCon.mat'));
    savepath = fullfile(folder_cell{i},save_folder);
    name = strcat(name_cell{i},'_Laguerre_DeCon.mat');
    
    % get gain information
    [gain_factor, gain_list] = gainCorrection(data);

    
    [~,non_decon_idx,full_data_size]= detectSaturation(data,low,high);
    %     handles.image_size(2) = full_data_size;
    %     fullIdx = 1:full_data_size(1);
    %     fullIdx=fullIdx';
    decon_idx = data.decon_idx;
    
%     if ~eq(num_of_batches,1)
           
    temp_text = sprintf('Deconvolving Line %d with %d points',i,length(decon_idx));
    set(handles.Decon_process_text,'String',temp_text);
    set(handles.hPb,'Value',(i)/num_of_batches*100);
    drawnow
    
%     end
    % reduced_data = out;
    % reduced_data(idx,:)=[];
    % currentdata  = get(data,'rawdata');
    % dif = currentdata-out;
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
    
    %     SNR1 = [];
    %     SNR2 = [];
    %     SNR3 = ones(full_data_size(1),1)*50;
    %     SNR4 = ones(full_data_size(1),1)*50;
    
    BG_scale_factor = zeros(full_data_size(1),1);
    % Tau1 = zeros(full_data_size(1),2);
    % Tau2 = zeros(full_data_size(1),2);
    CH1_LCs = zeros(full_data_size(1),12);
    CH2_LCs = zeros(full_data_size(1),12);
    CH3_LCs = zeros(full_data_size(1),12);
    CH4_LCs = zeros(full_data_size(1),12);
    % check whether there is data to be deconvolve
    if ~isempty(decon_idx)
        
        % curser_low = 1000;
        % curser_high = 1200;
        if ~isempty(BG_remove_postion)
        [BG_removed, BG_scale_factor] = removeBG(data,BG_remove_postion(1),BG_remove_postion(2));
        end
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
            SNR_temp = interp1(timeStamp,SNR_temp,Time);
            SNR{1} = SNR_temp;
            clear SNR_temp
            
            LaguerreModel_CH1=LaguerreModel(rawdatatruct1);
            LaguerreModel_CH1.iIRF_align();
            LaguerreModel_CH1.estimate_laguerre();
            handles.all_data{i,1}=LaguerreModel_CH1;
            
            LT1(decon_idx,:) = LaguerreModel_CH1.LTs';
            LT1 = interp1(timeStamp,LT1,Time);
            lifet_avg{1}=LT1;
            % Weight1(decon_idx,:) = expsLaguerreModel_CH1.Weights';
            INT1(decon_idx,:) = LaguerreModel_CH1.INTs';
            INT1 = INT1.*gain_factor;
            INT1 = interp1(timeStamp,INT1,Time);
            spec_int{1}=INT1;
            CH1_LCs(decon_idx,:) = LaguerreModel_CH1.LCs';
            Laguerre_coeffs{1} = CH1_LCs;
            alpha_out{1} = get(LaguerreModel_CH1,'alpha');
            
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
            SNR_temp = interp1(timeStamp,SNR_temp,Time);
            SNR{2} = SNR_temp;
            clear SNR_temp
            LaguerreModel_CH2=LaguerreModel(rawdatatruct2);
            LaguerreModel_CH2.iIRF_align();
            LaguerreModel_CH2.estimate_laguerre();
            handles.all_data{i,2}=LaguerreModel_CH2;
            
            LT2(decon_idx,:) = LaguerreModel_CH2.LTs';
            LT2 = interp1(timeStamp,LT2,Time);
            
            lifet_avg{2}=LT2;
            
            % Weight2(decon_idx,:) = expstruct2.Weights';
            INT2(decon_idx,:) = LaguerreModel_CH2.INTs';
            INT2 = INT2.*gain_factor;
            INT2 = interp1(timeStamp,INT2,Time);
            spec_int{2}=INT2;
            CH2_LCs(decon_idx,:) = LaguerreModel_CH2.LCs';
            Laguerre_coeffs{2} = CH2_LCs;
            alpha_out{2} = get(LaguerreModel_CH2,'alpha');
            
        end
        
        %% CH3 deconvolution
        
        if ~isempty(datastruct{3})
            rawdatatruct3 = channeldata(datastruct{3}',irf{3},dt,BW,decon_idx,noise,gain_list);
            SNR_temp = zeros(full_data_size(1),1);
            SNR_temp(decon_idx,:) = rawdatatruct3.SNR;
            SNR_temp = interp1(timeStamp,SNR_temp,Time);
            SNR{3} = SNR_temp;
            clear SNR_temp
            LaguerreModel_CH3=LaguerreModel(rawdatatruct3);
            LaguerreModel_CH3.iIRF_align();
            LaguerreModel_CH3.estimate_laguerre();
            handles.all_data{i,3}=LaguerreModel_CH3;
            
            LT3(decon_idx,:) = LaguerreModel_CH3.LTs';
            LT3 = interp1(timeStamp,LT3,Time);
            
            lifet_avg{3}=LT3;
            
            % Weight1(decon_idx,:) = expsLaguerreModel_CH1.Weights';
            INT3(decon_idx,:) = LaguerreModel_CH3.INTs';
            INT3 = INT3.*gain_factor;
            INT3 = interp1(timeStamp,INT3,Time);
            spec_int{3}=INT3;
            CH3_LCs(decon_idx,:) = LaguerreModel_CH3.LCs';
            Laguerre_coeffs{3} = CH3_LCs;
            alpha_out{3} = get(LaguerreModel_CH3,'alpha');
            
        end
        
        
        %% CH4 deconvolution
        if ~isempty(datastruct{4})
            rawdatatruct4 = channeldata(datastruct{4}',irf{4},dt,BW,decon_idx,noise,gain_list);
            SNR_temp = zeros(full_data_size(1),1);
            SNR_temp(decon_idx,:) = rawdatatruct4.SNR;
            SNR_temp = interp1(timeStamp,SNR_temp,Time);
            SNR{4} = SNR_temp;
            clear SNR_temp
            LaguerreModel_CH4=LaguerreModel(rawdatatruct4);
            LaguerreModel_CH4.iIRF_align();
            LaguerreModel_CH4.estimate_laguerre();
            handles.all_data{i,4}=LaguerreModel_CH4;
            
            LT4(decon_idx,:) = LaguerreModel_CH4.LTs';
            LT4 = interp1(timeStamp,LT4,Time);
            lifet_avg{4}=LT4;
            % Weight1(decon_idx,:) = expsLaguerreModel_CH1.Weights';
            INT4(decon_idx,:) = LaguerreModel_CH4.INTs';
            INT4 = INT4.*gain_factor;
            INT4 = interp1(timeStamp,INT4,Time);
            spec_int{4}=INT4;
            CH4_LCs(decon_idx,:) = LaguerreModel_CH4.LCs';
            Laguerre_coeffs{4} = CH4_LCs;
            alpha_out{4} = get(LaguerreModel_CH4,'alpha');
            
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
    
    
    INT_map{1} = cell2mat(INT1_map);
    INT_map{2} = cell2mat(INT2_map); 
    INT_map{3} = cell2mat(INT3_map);
    INT_map{4} = cell2mat(INT4_map);
    setappdata(handles.figure1,'INT_map',INT_map);
    
    
    LT_map{1} = cell2mat(LT1_map);
    LT_map{2} = cell2mat(LT2_map);
    LT_map{3} = cell2mat(LT3_map);
    LT_map{4} = cell2mat(LT4_map);
    setappdata(handles.figure1,'LT_map',LT_map);
    
    SNR_map{1} = cell2mat(SNR1_map);
    SNR_map{2} = cell2mat(SNR2_map);
    SNR_map{3} = cell2mat(SNR3_map);
    SNR_map{4} = cell2mat(SNR4_map);
    setappdata(handles.figure1,'SNR_map',SNR_map);
    
    % save to .mat file
    cd(savepath)
    save(name,'lifet_avg','spec_int','alpha_out','Laguerre_coeffs',...
        'timeStamp','SNR','test','decon_idx','full_data_size','BG',...
        'BG_scale_factor','gain_factor','gain_list','Data');
%     
    
    
end
    cd(handles.rootFolder)