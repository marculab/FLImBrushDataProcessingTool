function [data,gain_values,timeStamp,spaceStamp,pos_values,line_indx,ID]= LoadRawTRIPLEXData(nameListIn,folderListIn)
% nameListIn: 1D cell array of filename
% folderListIn: 1D cell array of filepath

folderListIn = folderListIn{1};

if ~isempty(nameListIn)
    % check total data size
    nfiles = length(nameListIn); %number of data files
    % filename_temp =cell2mat(nameListIn(1));
    % [~, data_temp] = DataLoader(fullfile(folderListIn, filename_temp));
    % temp = whos('data_temp');
    % filesize = temp.bytes*1e-6; % single file size in MB
    % total_data_size = filesize*num_of_files;
    % Message = sprintf('The estimated total data size for deconvolution is %.2f MB',total_data_size);
    % h = msgbox(Message);
    % waitfor(h);
    
    clear data_temp header_temp temp filesize total_data_size Message
    % initilize data matrix
    filenames = cell(nfiles, 1);
    prefixes = cell(nfiles, 1);
    % gain_values = cell(num_of_files, 1);
    % timeStamp = cell(num_of_files, 1);
    % pos_values = cell(num_of_files, 1);
    line_looking_table = zeros(nfiles, 1); % stores number of measurement 
    % per data file, used to recover coresponding data points later.
    
    % initilize all variables to empty
    data = [];
%     header = [];
    gain_values = [];
    timeStamp = [];
    pos_values = [];
    spaceStamp = [];
    ID = [];
    for iFiles=1:nfiles
        nameStr = cell2mat(nameListIn(iFiles));
        line_num = str2double(nameStr(end-2:end));
        filenames{iFiles, 1} = fullfile(folderListIn, nameStr);
        prefixes{iFiles, 1} = nameStr;
        
        % load in both header and raw data
        [header_temp, data_temp] = DataLoader(fullfile(folderListIn, nameStr)); % load data
        %remove noise(DC)
        % calculate noise from all points last 100 data point
        %     noise = mean(data_temp(:,(end - 100):end),2);
        % calculate noise from all points first 100 data point
        %     noise = mean(data_temp(:,1:100),2);
        %     data_temp = bsxfun(@minus, data_temp, noise);
        
        number_of_points_temp = size(data_temp,2); % get number of measurement point in this file
        
        line_looking_table(iFiles) = number_of_points_temp; % store number of measurement points
        
        try
            gain_values_temp = header_temp.GainVoltage; % check whether this is gain information
        catch err
            gain_values_temp = zeros(number_of_points_temp,1); % if no gain information, set all to zero
        end
        
        try
            timeStamp_temp = header_temp.TimeStamp; % check time stamp information
        catch err
            timeStamp_temp = zeros(number_of_points_temp,1); % if no time stamp, set all to zero
        end
        
        try
            pos_values_temp = header_temp.FiberPos; % check location information for hand scanning system
        catch err
            pos_values_temp = zeros(number_of_points_temp,1); % if no information available, set all to zero
        end
        
        try
            spaceStamp_temp = header_temp.SpaceStamp; % check spacestamp info, will be deprated in next release
        catch err
            spaceStamp_temp = zeros(number_of_points_temp,1); % set spacestamp to zeros
        end
        
        try ID_temp = header_temp.ID; % check ID info
        catch err
            ID_temp = []; % if no ID info, set to empty
        end
        
        %fix empty space stamp
        %     if all(~header_temp.SpaceStamp)||(nnz(header_temp.SpaceStamp)<5)
        %         try
        %         spaceStamp_temp = load(fullfile(folderListIn, [str '_SpaceStamp.txt']));
        %         catch err
        %
        %         end
        %     end
        
        data = [data data_temp]; % append new data to end
        gain_values = [gain_values; gain_values_temp]; % append gain 
        timeStamp = [timeStamp; timeStamp_temp]; % append timeStamp
        spaceStamp = [spaceStamp; spaceStamp_temp]; % append spaceStamo
        pos_values = [pos_values; pos_values_temp]; % appen position info
        
        line_indx(iFiles,1) = sum(line_looking_table(1:iFiles)); % record starding line index of each file
        
        ID = [ID ID_temp]; % append ID
        
    end
    % set data to positive
    data = -data;
    
    if isfield(header_temp,'SScanFlag')
        if header_temp.SScanFlag % check if it is a s scan
            if ~mod(line_num,2) % check line number, if it is even, flip
                data = flipud(data);
                gain_values = flipud(gain_values);
                timeStamp = flipud(timeStamp);
                spaceStamp = flipud(spaceStamp);
                pos_values = flipud(pos_values);
            end
        end
    end
%     return
else
%     return
end

end

