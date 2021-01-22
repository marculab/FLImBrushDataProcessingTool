function handles = run_Decon(Decon_param,Files,Background,iIRF,handles)
namelist = Files.namelist(cell2mat(Files.selector);
folderlist = Files.folderlist(cell2mat(Files.selector);
num_of_batches = length(namelist);

t=datetime('now');
formatOut = 'mmm-dd-yyyy HH.MM PM';

% so far only LMEA method go through this processor
save_folder = ['DeCon_LaguerreMEA ' datestr(t,formatOut)];
mkdir(folderlist{1},save_folder); % make folder to save data


raw_data_struct_cell = cell(size(namelist));
% resultMaps = struct('INT_maps', cell(4,1), ...
%     'INT_maps_interped', cell(4,1), ...
%     'LT_maps', cell(4,1), ...
%     'LT_maps_interped', cell(4,1), ...
%     'SNR_maps', cell(4,1), ...
%     'SNR_maps_interped', cell(4,1), ...
%     'Laguerre_maps', cell(4,1), ...
%     'Laguerre_maps_interped', cell(4,1), ...
%     'exp_maps', cell(4,1), ...
%     'exp_maps_interped', cell(4,1), ...
%     'LMEA_maps', cell(4,1), ...
%     'LMEA_maps', cell(4,1) ...
%     );

INT1_map  = cell(size(namelist));
LT1_map  = cell(size(namelist));
INT2_map  = cell(size(namelist));
LT2_map  = cell(size(namelist));
INT3_map  = cell(size(namelist));
LT3_map  = cell(size(namelist));
INT4_map  = cell(size(namelist));
LT4_map  = cell(size(namelist));

INT1_map_interped  = cell(size(namelist));
LT1_map_interped  = cell(size(namelist));
INT2_map_interped  = cell(size(namelist));
LT2_map_interped  = cell(size(namelist));
INT3_map_interped  = cell(size(namelist));
LT3_map_interped  = cell(size(namelist));
INT4_map_interped  = cell(size(namelist));
LT4_map_interped  = cell(size(namelist));

SNR1_map  = cell(size(namelist));
SNR2_map  = cell(size(namelist));
SNR3_map  = cell(size(namelist));
SNR4_map  = cell(size(namelist));

SNR1_map_interped  = cell(size(namelist));
SNR2_map_interped  = cell(size(namelist));
SNR3_map_interped  = cell(size(namelist));
SNR4_map_interped  = cell(size(namelist));

CH1_Taus_1_map = cell(size(namelist));
CH2_Taus_1_map = cell(size(namelist));
CH3_Taus_1_map = cell(size(namelist));
CH4_Taus_1_map = cell(size(namelist));

CH1_Taus_2_map = cell(size(namelist));
CH2_Taus_2_map = cell(size(namelist));
CH3_Taus_2_map = cell(size(namelist));
CH4_Taus_2_map = cell(size(namelist));

CH1_Taus_1_map_interped = cell(size(namelist));
CH2_Taus_1_map_interped = cell(size(namelist));
CH3_Taus_1_map_interped = cell(size(namelist));
CH4_Taus_1_map_interped = cell(size(namelist));

CH1_Taus_2_map_interped = cell(size(namelist));
CH2_Taus_2_map_interped = cell(size(namelist));
CH3_Taus_2_map_interped = cell(size(namelist));
CH4_Taus_2_map_interped = cell(size(namelist));

CH1_Weight_1_map = cell(size(namelist));
CH2_Weight_1_map = cell(size(namelist));
CH3_Weight_1_map = cell(size(namelist));
CH4_Weight_1_map = cell(size(namelist));

CH1_Weight_2_map = cell(size(namelist));
CH2_Weight_2_map = cell(size(namelist));
CH3_Weight_2_map = cell(size(namelist));
CH4_Weight_2_map = cell(size(namelist));

CH1_Weight_1_map_interped = cell(size(namelist));
CH2_Weight_1_map_interped = cell(size(namelist));
CH3_Weight_1_map_interped = cell(size(namelist));
CH4_Weight_1_map_interped = cell(size(namelist));

CH1_Weight_2_map_interped = cell(size(namelist));
CH2_Weight_2_map_interped = cell(size(namelist));
CH3_Weight_2_map_interped = cell(size(namelist));
CH4_Weight_2_map_interped = cell(size(namelist));

% main GUI variable to save all the data structure
all_data = cell(length(namelist),4);

% set progress bar to 0
set(handles.hPb,'Value',0);

for i=1:num_of_batches

    %% initilize all local variables
    % initilize common variables
    lifet_avg = cell(4,1);  % average lifetime
    lifet_avg_interped = cell(4,1); % interpreted average lifetime
    spec_int = cell(4,1);   % integrated intensity
    spec_int_interped = cell(4,1);   % integrated intensity
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
    alpha_out = cell(4,1);  % alpha value for each channel
    Laguerre_coeffs = cell(4,1);    % Laguerre coefficinets
    Weights_LMEA = cell(4,1);
    %% preprocessing data
    
    rawDataClass = RawDataStrcutClass(namelist(i),folderlist(i),Background.data,iIRF.data); % creat raw data structure class
    n_pix = rawDataClass.rawdata_dimention(1); % get number of pixels per line
    ID = get(rawDataClass,'ID'); % ID of measurement
    timeStamp = get(rawDataClass,'timeStamp'); % get time stamp for interpoltion of the data
    spaceStamp = get(rawDataClass,'spaceStamp');% get space stmap for interpolation of data
    %     DAQtime = round(mean(diff(timeStamp))*1000); % calculate average DAQ time
    if all(~spaceStamp)
    Time = (0:n_pix-1)*Decon_param.DAQ_time/1000; % acq time points for a frame
    Time = Time'; % reshape time
    else
        Time = spaceStamp;
    end
    peak_position = iIRF.peaks.*Decon_param.Channels;

    [pathstr,name,ext] = fileparts(namelist{i}); % get file info
    savepath = fullfile(folderlist{i},save_folder); % set saving pass
    name = strcat(namelist{i},'_Laguerre_DeCon.mat'); % creat saving name
    
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
    BG_scale_factor = zeros(full_data_size(1),1);
    CH_LCs = zeros(full_data_size(1),12);

    % check whether there is data to be deconvolve
    if ~isempty(decon_idx)

        % curser_low = 1000;
        % curser_high = 1200;


        if ~isempty(Background.refRange)
            [BG_removed, BG_scale_factor] = removeBG(rawDataClass,Background.refRange(1),Background.refRange(2));

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
                rawdatastruct = channeldata(datastruct{j}',iIRF.data{j},Decon_param.time_res,Decon_param.bw,decon_idx,noise,gain_list);
                SNR_temp = zeros(full_data_size(1),1);
                SNR_temp(decon_idx,:) = rawdatastruct.SNR;
                SNR{j} = SNR_temp;
                SNR_temp = interp1(timeStamp,SNR_temp,Time);
                SNR_interped{j} = SNR_temp;
                LT = zeros(full_data_size(1),1);
                INT = zeros(full_data_size(1),1);
                
                Laguerre_Struct=LaguerreModel(rawdatastruct);
                Laguerre_Struct.iIRF_align();
                Laguerre_Struct.estimate_laguerre();
                Data{j} = Laguerre_Struct;
                all_data{i,j}=Laguerre_Struct;
                
                LT(decon_idx,:) = Laguerre_Struct.LTs';
                LT_interped = interp1(timeStamp,LT,Time);
                lifet_avg{j}=LT;
                lifet_avg_interped{j}=LT_interped;
                
                % Weight1(decon_idx,:) = expsLaguerreModel_CH1.Weights';
                INT(decon_idx,:) = Laguerre_Struct.INTs';
                INT = INT.*gain_factor;
                INT_interped = interp1(timeStamp,INT,Time);
                spec_int{j}=INT;
                spec_int_interped{j}=INT_interped;
                
                CH_LCs(decon_idx,:) = Laguerre_Struct.LCs';
                Laguerre_coeffs{j} = CH_LCs;
                alpha_out{j} = get(Laguerre_Struct,'alpha');
                
                % clear variables for each line
                clear SNR_temp LT LT_interped INT INT_interped CH_LCs LaguerreStruct Taus Taus_interped Weights Weights_interped
            end

        end



    else % if nothing to deconvolve, set everything to 0

        for k=1:4
            spec_int{k} = zeros(full_data_size(1),1);
            lifet_avg{k} = zeros(full_data_size(1),1);
            SNR{k} = zeros(full_data_size(1),1);
            spec_int_interped{k} = zeros(full_data_size(1),1);
            lifet_avg_interped{k} = zeros(full_data_size(1),1);
            SNR_interped{k} = zeros(full_data_size(1),1);
        end
    end
    %% save data to cell arrays
    INT1_map{i}  = spec_int{1}';
    LT1_map{i}  = lifet_avg{1}';
    SNR1_map{i}  = SNR{1}';
    INT2_map{i} = spec_int{2}';
    LT2_map{i}  = lifet_avg{2}';
    SNR2_map{i}  = SNR{2}';
    INT3_map{i}  = spec_int{3}';
    LT3_map{i}  = lifet_avg{3}';
    SNR3_map{i}  = SNR{3}';
    INT4_map{i} = spec_int{4}';
    LT4_map{i}  = lifet_avg{4}';
    SNR4_map{i}  = SNR{4}';

    INT1_map_interped{i}  = spec_int_interped{1}';
    LT1_map_interped{i}  = lifet_avg_interped{1}';
    SNR1_map_interped{i}  = SNR_interped{1}';
    INT2_map_interped{i} = spec_int_interped{2}';
    LT2_map_interped{i}  = lifet_avg_interped{2}';
    SNR2_map_interped{i}  = SNR_interped{2}';
    INT3_map_interped{i}  = spec_int_interped{3}';
    LT3_map_interped{i}  = lifet_avg_interped{3}';
    SNR3_map_interped{i}  = SNR_interped{3}';
    INT4_map_interped{i} = spec_int_interped{4}';
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
    cd(savepath)
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

end % end of the loop of line


end