function handles = exportPreprossedData(Decon_method, name_cell, folder_cell,BG,irf,handles,show_waveform_flag)

low = getappdata(handles.figure1,'low');
high = getappdata(handles.figure1,'high');
dt = getappdata(handles.figure1,'dt');
BW = getappdata(handles.figure1,'BW');
truncation_length = getappdata(handles.figure1,'truncation_length');
num_of_batches = length(name_cell);
save_folder = 'PreprocessedData';
mkdir(folder_cell{1},save_folder);
rawDataAll = [];
truncatedDataAll = [];
CH1DataAll=[];
CH2DataAll=[];
CH3DataAll=[];
CH4DataAll=[];

%% loop for each line
for i=1:num_of_batches
    
    
    
    %% preprocessing data
    data = RawDataStrcutClass(name_cell(i),folder_cell(i),BG,irf); % creat raw data structure class
    n_pix = data.rawdata_dimention(1); % get number of pixels per line
    ID = get(data,'ID'); % ID of measurement
    timeStamp = get(data,'timeStamp'); % get time stamp for interpoltion of the data
    spaceStamp = get(data,'spaceStamp');% get space stmap for interpolation of data
    DAQtime = str2double(get(handles.DAQ_time,'String'));
    %     DAQtime = round(mean(diff(timeStamp))*1000); % calculate average DAQ time
    if all(~spaceStamp)
        Time = (0:n_pix-1)*DAQtime/1000; % acq time points for a frame
        Time = Time'; % reshape time
    else
        Time = spaceStamp;
    end
    BG_remove_postion = getappdata(handles.figure1,'BG_remove_postion'); % ger BG subtraction position
    peak_position = getappdata(handles.figure1,'peak_position'); % get deconvolution peak position
    channel_selection = getappdata(handles.figure1,'channel_selection'); % get deconvolution channel_selection
    peak_position = peak_position.*channel_selection;
    
    [pathstr,name,ext] = fileparts(name_cell{i}); % get file info
    savepath = fullfile(folder_cell{i},save_folder); % set saving pass
    
    % get gain information
    [gain_factor, gain_list] = gainCorrection(data);
    
    
    noise = get_noise(data); %get noise form data
    %% initialize data matrix
    
    
    if ~isempty(BG_remove_postion)
        [~, BG_scale_factor] = removeBG(data,BG_remove_postion(1),BG_remove_postion(2));
    end
    
    truncatedData = truncateData(data,truncation_length,peak_position); %cell array to store truncated data
    CH1DataAll = cat(3,CH1DataAll,truncatedData{1});
    CH2DataAll = cat(3,CH2DataAll,truncatedData{2});
    CH3DataAll = cat(3,CH3DataAll,truncatedData{3});
    CH4DataAll = cat(3,CH4DataAll,truncatedData{4});
    
    truncatedDataArray = cell2mat(truncatedData);
    truncatedDataAll = cat(3,truncatedDataAll,truncatedDataArray);
    rawData  = get(data,'rawdata');
    BG_removed = get(data,'data');
    
    rawDataAll = [rawDataAll;BG_removed];
    
    
    
    % save to .mat file
    cd(savepath)
    
    name = strcat(name_cell{i},'Preprocessed');
    
    save(name,'rawData','BG_removed','truncatedData','timeStamp','spaceStamp','gain_factor', ...
        'gain_list','truncatedDataArray');
    
end
save(strcat(name,'All'),'truncatedDataAll','CH1DataAll','CH2DataAll','CH3DataAll','CH4DataAll');
h = msgbox('Finished output');

cd(handles.rootFolder)