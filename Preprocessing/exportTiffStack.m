function handles = exportTiffStack(Decon_method, name_cell, folder_cell,BG,irf,handles,show_waveform_flag)

low = getappdata(handles.figure1,'low');
high = getappdata(handles.figure1,'high');
dt = getappdata(handles.figure1,'dt');
BW = getappdata(handles.figure1,'BW');
truncation_length = getappdata(handles.figure1,'truncation_length');
num_of_batches = length(name_cell);
save_folder = 'tiffStack';
mkdir(folder_cell{1},save_folder);
rawDataAll = [];
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
    
    %     [~,non_decon_idx,full_data_size]= detectSaturation(data,low,high);
    
    %     decon_idx = data.decon_idx; % get deconvolution index
    noise = get_noise(data); %get noise form data
    %% initialize data matrix
    
    
    if ~isempty(BG_remove_postion)
        [~, BG_scale_factor] = removeBG(data,BG_remove_postion(1),BG_remove_postion(2));
    end
    
    
    % handles.truncation_length=680;
    % peak_position =[515 1281 0 0];
    
    truncatedData = truncateData(data,truncation_length,peak_position); %cell array to store truncated data
    truncatedDataArray = cell2mat(truncatedData);
%     rawData  = get(data,'rawdata');
%     BG_removed = get(data,'data');
    
    rawDataAll = [rawDataAll;truncatedDataArray];
    
    
    
    %% save to .mat file
%     cd(savepath)
    
%     name = strcat(name_cell{i},'Preprocessed');
    
%     save(name,'rawData','BG_removed','truncatedData','timeStamp','spaceStamp','gain_factor', 'gain_list');
    
end
cd(savepath)
[r,c] = size(rawDataAll);
nlay  = num_of_batches;
rawdataCube   = permute(reshape(rawDataAll',[c,r/nlay,nlay]),[2,3,1]);
% save('ZStack.mat','rawdataCube');

imwrite(rawdataCube(:,:,1), 'imgstack.tif')
for k = 2:size(rawdataCube,3)
    imwrite(rawdataCube(:,:,k), 'imgstack.tif','tif', 'WriteMode', 'append');
end

% imwrite(rawDataAll(:,1), 'rowstack.tif')
% for k = 2:size(rawDataAll,2)
%     imwrite(rawDataAll(:,k), 'rowstack.tif','tif', 'WriteMode', 'append');
% end
h = msgbox('Finished output');
cd(handles.rootFolder)