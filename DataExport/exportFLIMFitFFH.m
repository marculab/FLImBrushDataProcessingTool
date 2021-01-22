function handles = exportFLIMFitFFH(Decon_param, name_cell, folder_cell,Background,iIRF,handles,show_waveform_flag)
%% get necessary parameter
num_of_batches = length(name_cell);

%% get current date and time and generate saving folder name
t=datetime('now');
formatOut = 'mmm-dd-yyyy HH.MM PM';


save_folder = ['FLIMFit Compatible ffh ' datestr(t,formatOut)];


mkdir(folder_cell{1},save_folder); % make folder to save data

imgName = name_cell{1};
imgName(end-2:end) = [];
savename = [imgName '.ffh'];
%% set progress bar to 0
set(handles.hPb,'Value',0);

%%initilize data matrix
BGRemovedRawDataAll = [];
truncatedDataAll = [];
CH1DataAll=[];
CH2DataAll=[];
CH3DataAll=[];
CH4DataAll=[];
timeStampImage = [];
DAQ_time = Decon_param.DAQ_time;

%% loop for each line
for i=1:num_of_batches
    %% preprocessing data
    
    rawDataClass = RawDataStrcutClass(name_cell(i),folder_cell(i),Background.data,iIRF.data); % creat raw data structure class
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
    
    [pathstr,name,ext] = fileparts(name_cell{i}); % get file info
    savepath = fullfile(folder_cell{i},save_folder); % set saving pass
    
%     name = strcat(name_cell{i},'Preprocessed Data'); % creat saving file name
    
    % get gain information
    
    [gain_factor, gain_list] = gainCorrection(rawDataClass);
    
    noise = get_noise(rawDataClass); %get noise form data
    
    %% set GUI progress bar
    temp_text = sprintf('Loading Line %d ',i);
    set(handles.edit_decon_process,'String',temp_text);
    set(handles.hPb,'Value',(i)/num_of_batches*100);
    drawnow
   
    if ~isempty(Background.refRange)
        [BG_removed, BG_scale_factor] = removeBG(rawDataClass,Background.refRange(1),Background.refRange(2));
        
    end
    
    if show_waveform_flag
        plot(BG_removed');
        ylim([-0.1 0.75])
        drawnow
    end
    
    
    truncateDataAllCell = truncateData(rawDataClass,Decon_param.ch_width,peak_position); %cell array to store truncated data
    
    CH1DataAll = cat(3,CH1DataAll,truncateDataAllCell{1});
    CH2DataAll = cat(3,CH2DataAll,truncateDataAllCell{2});
    CH3DataAll = cat(3,CH3DataAll,truncateDataAllCell{3});
    CH4DataAll = cat(3,CH4DataAll,truncateDataAllCell{4});
    
    timeStampImage = [timeStampImage;timeStamp];
    
    
    
    truncatedDataArray = cell2mat(truncateDataAllCell);
    truncatedDataAll = cat(3,truncatedDataAll,truncatedDataArray); % truncate data concascated
    
    rawData  = get(rawDataClass,'rawdata');
    
    BG_removed = get(rawDataClass,'data');
    BGRemovedRawDataAll = [BGRemovedRawDataAll;BG_removed]; % not truncated data
    
    
    %% save to .mat file
%     cd(savepath)
%     BG = Background.data;
    
%     save(name,'rawData','BG_removed','truncateDataAllCell','timeStamp','spaceStamp','gain_factor', ...
%         'gain_list','truncatedDataArray');
   
    clear datastruct rawDataClass timeStamp rawData BG_removed truncatedDataArray
end % end of the loop of line

%% save 3D datacube
DataCube(:,1,:,:) = CH1DataAll;
DataCube(:,2,:,:) = CH2DataAll;
DataCube(:,3,:,:) = CH3DataAll;
DataCube(:,4,:,:) = CH4DataAll;
DataCube = permute(DataCube,[3,2,4,1]);
% DataCube = cat(4,DataCube,DataCube);
DataCube = uint16(DataCube*1000);
t = (1:size(DataCube,1))*80;
cd(savepath)
writeFfh(savename, DataCube, t);

% save(strcat(name,'ImageCube'),'truncatedDataAll','CH1DataAll','CH2DataAll','CH3DataAll','CH4DataAll');
h = msgbox('Finished output');


cd(handles.rootFolder)
