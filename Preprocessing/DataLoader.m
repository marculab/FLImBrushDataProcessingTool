function [Header, Data] = DataLoader(filename)
% DATALOADER load data saved by Labview for multiple systems
%   [Header, Data] = DataLoader(path,filename)
%   path - string of the data/header directory.
%   filename - string of data file name, the header file name must be
%       path + filename + '_header'
%   Header - structure containing DAQ parameters, and other information
%       about the data set, e.g. time stamp, gain voltage, etc.
%   Data - Data matrix in "double" word width that corresponds to the
%       voltage signal at the digitizer input. 1st dimentsion is time.
%   This function requires dependency of:
%       function data=loadBinary(filepath,type,saveflag)
%
%   Create by Michael Ma on Nov 2014
%   last updated by Xiangnan on Apr 23rd 2019
%
Hpath = [filename '_header'];
Dpath = filename;
% SpaceStampPath = [filename '_SpaceStamp.txt'];
% TimeStampPath = [filename '_TimeStamp.txt'];
fh = fopen(Hpath);
if fh == -1
    error('DataLoader: cannot open header file!');
else
    SysID = fread(fh,[1 1],'int32');
    switch SysID
        % msTRFS Davis
        case 1
            InputRange = fread(fh,[1 1],'double');
            Offset = fread(fh,[1 1],'double');
            TimeResolution = fread(fh,[1 1],'double');
            BW = fread(fh,[1 1],'double');
            numAvg = fread(fh,[1 1],'int32');
            n = fread(fh,[1 1],'int32');
            m = fread(fh,[1 1],'int32');
            ts = fread(fh,[m,n],'double');
            TimeStamp = ts(:,1);
            GainVoltage = ts(:,2);
            Header = struct('SysID',SysID,'InputRange',InputRange,'Offset',Offset,...
                'TimeResolution',TimeResolution,'BW',BW,'numAvg',numAvg,...
                'TimeStamp',TimeStamp,'GainVoltage',GainVoltage);
            % msTRFS Sacramento V4
        case 2
            InputRange = fread(fh,[1 1],'double');
            Offset = fread(fh,[1 1],'double');
            TimeResolution = fread(fh,[1 1],'double');
            BW = fread(fh,[1 1],'double');
            numAvg = fread(fh,[1 1],'int32');
            n = fread(fh,[1 1],'int32');
            m = fread(fh,[1 1],'int32');
            ts = fread(fh,[m,n],'double');
            TimeStamp = ts(:,1);
            GainVoltage = ts(:,2);
            try
            ID = ts(:,3);
            catch
               ID = zeros(size(GainVoltage));
            end
%             ID = zeros(length(GainVoltage),1);
            Header = struct('SysID',SysID,'InputRange',InputRange,'Offset',Offset,...
                'TimeResolution',TimeResolution,'BW',BW,'numAvg',numAvg,...
                'TimeStamp',TimeStamp,'GainVoltage',GainVoltage,'ID',ID);
            % TRFS Sacramento
        case 3
            InputRange = fread(fh,[1 1],'double');
            Offset = fread(fh,[1 1],'double');
            TimeResolution = fread(fh,[1 1],'double');
            BW = fread(fh,[1 1],'double');
            numAvg = fread(fh,[1 1],'int32');
            n = fread(fh,[1 1],'int32');
            m = fread(fh,[1 1],'int32');
            ts = fread(fh,[m,n],'double');
            Wavelength = ts(:,1);
            GainVoltage = ts(:,2);
            Header = struct('SysID',SysID,'InputRange',InputRange,'Offset',Offset,...
                'TimeResolution',TimeResolution,'BW',BW,'numAvg',numAvg,...
                'Wavelength',Wavelength,'GainVoltage',GainVoltage);
            % FLIm catheter
        case 4
            InputRange = fread(fh,[1 1],'double');
            Offset = fread(fh,[1 1],'double');
            TimeResolution = fread(fh,[1 1],'double');
            BW = fread(fh,[1 1],'double');
            numAvg = fread(fh,[1 1],'int32');
            PBdistance = fread(fh,[1 1],'double');
            PBspeed = fread(fh,[1 1],'double');
            Gain = fread(fh,[1 1],'double');
            Catheter = fread(fh,[1 1],'uint8');
            n = fread(fh,[1 1],'int32');
            m = fread(fh,[1 1],'int32');
            TimeStamp = fread(fh,[m,n],'double');
            Header = struct('SysID',SysID,'InputRange',InputRange,'Offset',Offset,...
                'TimeResolution',TimeResolution,'BW',BW,'numAvg',numAvg,...
                'Gain',Gain,'Catheter',Catheter,...
                'PBdistance',PBdistance,'PBspeed',PBspeed,'TimeStamp',TimeStamp);
            %         'Gain',Gain,'Catheter',Catheter,...
            % IVUS catheter
        case 5
            InputRange = fread(fh,[1 1],'double');
            Offset = fread(fh,[1 1],'double');
            TimeResolution = fread(fh,[1 1],'double');
            n = fread(fh,[1 1],'int32');
            m = fread(fh,[1 1],'int32');
            TimeStamp = fread(fh,[m,n],'double');
            TimeStamp = TimeStamp(1:floor(length(TimeStamp)/256)*256);
            Header = struct('SysID',SysID,'InputRange',InputRange,'Offset',Offset,...
                'TimeResolution',TimeResolution,'TimeStamp',TimeStamp);
            % Davis TRFS
        case 6
            InputRange = fread(fh,[1 1],'double');
            Offset = fread(fh,[1 1],'double');
            TimeResolution = fread(fh,[1 1],'double');
            BW = fread(fh,[1 1],'double');
            numAvg = fread(fh,[1 1],'int32');
            n = fread(fh,[1 1],'int32');
            m = fread(fh,[1 1],'int32');
            ts = fread(fh,[m,n],'double');
            Wavelength = ts(:,1);
            GainVoltage = ts(:,2);
            Header = struct('SysID',SysID,'InputRange',InputRange,'Offset',Offset,...
                'TimeResolution',TimeResolution,'BW',BW,'numAvg',numAvg,...
                'Wavelength',Wavelength,'GainVoltage',GainVoltage);
            % Davis msTRFS with fiber motion
        case 9
            InputRange = fread(fh,[1 1],'double');
            Offset = fread(fh,[1 1],'double');
            TimeResolution = fread(fh,[1 1],'double');
            BW = fread(fh,[1 1],'double');
            numAvg = fread(fh,[1 1],'int32');
            n = fread(fh,[1 1],'int32');
            m = fread(fh,[1 1],'int32');
            ts = fread(fh,[m,n],'double');
            TimeStamp = ts(:,1);
            GainVoltage = ts(:,2);
            Position = ts(:,3);
            Header = struct('SysID',SysID,'InputRange',InputRange,'Offset',Offset,...
                'TimeResolution',TimeResolution,'BW',BW,'numAvg',numAvg,...
                'TimeStamp',TimeStamp,'GainVoltage',GainVoltage,'Position',Position);
            % TRFS Davis NI-Horiba
        case 26
            InputRange = fread(fh,[1 1],'double');
            Offset = fread(fh,[1 1],'double');
            TimeResolution = fread(fh,[1 1],'double');
            BW = fread(fh,[1 1],'double');
            numAvg = fread(fh,[1 1],'int32');
            n = fread(fh,[1 1],'int32');
            m = fread(fh,[1 1],'int32');
            ts = fread(fh,[m,n],'double');
            Wavelength = ts(:,1);
            GainVoltage = ts(:,2);
            Header = struct('SysID',SysID,'InputRange',InputRange,'Offset',Offset,...
                'TimeResolution',TimeResolution,'BW',BW,'numAvg',numAvg,...
                'Wavelength',Wavelength,'GainVoltage',GainVoltage);
            % Sacramento TRFS with Chromex Monochromator/Old software
        case {31,32}
            InputRange = fread(fh,[1 1],'double');
            Offset = fread(fh,[1 1],'double');
            TimeResolution = fread(fh,[1 1],'double');
            BW = fread(fh,[1 1],'double');
            numAvg = fread(fh,[1 1],'int32');
            n = fread(fh,[1 1],'int32');
            m = fread(fh,[1 1],'int32');
            ts = fread(fh,[m,n],'double');
            Wavelength = ts(:,1);
            GainVoltage = ts(:,2);
            Header = struct('SysID',SysID,'InputRange',InputRange,'Offset',Offset,...
                'TimeResolution',TimeResolution,'BW',BW,'numAvg',numAvg,...
                'Wavelength',Wavelength,'GainVoltage',GainVoltage);
            % MMC legacy format
        case 33
            InputRange = fread(fh,[1 1],'double');
            Offset = fread(fh,[1 1],'double');
            TimeResolution = fread(fh,[1 1],'double');
            BW = fread(fh,[1 1],'double');
            numAvg = fread(fh,[1 1],'int32');
            PBdistance = fread(fh,[1 1],'double');
            PBspeed = fread(fh,[1 1],'double');
            Gain = fread(fh,[1 1],'double');
            Catheter = fread(fh,[1 1],'uint8');
            Header = struct('SysID',SysID,'InputRange',InputRange,...
                'Offset',Offset,'TimeResolution',TimeResolution,...
                'BW',BW,'numAvg',numAvg,'Gain',Gain,...
                'Catheter',Catheter,'PBdistance',PBdistance,...
                'PBspeed',PBspeed);
            % FLIm scan legacy format
        case 34
            InputRange = fread(fh,[1 1],'double');
            Offset = fread(fh,[1 1],'double');
            TimeResolution = fread(fh,[1 1],'double');
            Header = struct('SysID',SysID,'InputRange',InputRange,'Offset',Offset,...
                'TimeResolution',TimeResolution);
        case 35 %TE Data V5, raster scan
            InputRange = fread(fh,[1 1],'double');
            Offset = fread(fh,[1 1],'double');
            TimeResolution = fread(fh,[1 1],'double');
            BW = fread(fh,[1 1],'double');
            numAvg = fread(fh,[1 1],'int32');
            n = fread(fh,[1 1],'int32');
            m = fread(fh,[1 1],'int32');
            ts = fread(fh,[m,n],'double');
            TimeStamp = ts(:,1);
            GainVoltage = ts(:,2);
            SScanFlag = 0;% flag for S scan
            try
                SpaceStamp = ts(:,3);
            catch
                SpaceStamp = zeros(size(TimeStamp));
            end
            Header = struct('SysID',SysID,'InputRange',InputRange,'Offset',Offset,...
                'TimeResolution',TimeResolution,'BW',BW,'numAvg',numAvg,...
                'TimeStamp',TimeStamp,'GainVoltage',GainVoltage,'SpaceStamp',SpaceStamp,'SScanFlag',SScanFlag);
        case 36 %TE Data V5 with space stamp
            InputRange = fread(fh,[1 1],'double');
            Offset = fread(fh,[1 1],'double');
            TimeResolution = fread(fh,[1 1],'double');
            BW = fread(fh,[1 1],'double');
            numAvg = fread(fh,[1 1],'int32');
            n = fread(fh,[1 1],'int32');
            m = fread(fh,[1 1],'int32');
            ts = fread(fh,[m,n],'double');
            TimeStamp = ts(:,1);
            GainVoltage = ts(:,2);
            SScanFlag = 1;% flag for S scan
            try
                SpaceStamp = ts(:,3);
            catch
                SpaceStamp = zeros(size(TimeStamp));
            end
            Header = struct('SysID',SysID,'InputRange',InputRange,'Offset',Offset,...
                'TimeResolution',TimeResolution,'BW',BW,'numAvg',numAvg,...
                'TimeStamp',TimeStamp,'GainVoltage',GainVoltage,'SpaceStamp',SpaceStamp,'SScanFlag',SScanFlag);
        case 101
            InputRange = fread(fh,[1 1],'double');
            Offset = fread(fh,[1 1],'double');
            TimeResolution = fread(fh,[1 1],'double');
            BW = fread(fh,[1 1],'double');
            numAvg = fread(fh,[1 1],'int32');
            n = fread(fh,[1 1],'int32');
            m = fread(fh,[1 1],'int32');
            ts = fread(fh,[m,n],'double');
            TimeStamp = ts(:,1);
            GainVoltage = ts(:,2);
            Header = struct('SysID',SysID,'InputRange',InputRange,'Offset',Offset,...
                'TimeResolution',TimeResolution,'BW',BW,'numAvg',numAvg,...
                'TimeStamp',TimeStamp,'GainVoltage',GainVoltage);
            % unknown
        case default
            error('DataLoader: unknown system ID! No header will be loaded.');
    end
    fclose(fh);
    fd = fopen(Dpath);
    if fd == -1
        error('DataLoader: cannot open data file!');
    else
        fclose(fd);
        % load Data
        switch SysID
            case 1
                Data = loadBinary(Dpath,'int16',0);
                Data = Data/256/numAvg*InputRange + Offset;
            case 2
                [~,Data] = loadBinary(Dpath,'double',0); % not all data is
                % retrived from Raw data files. Only those processed real
                % time is retrived. Raw data has more data point than
                % processed data
            case 3
                Data = loadBinary(Dpath,'double',0);
            case 4
                Data = loadBinary(Dpath,'int16',0);
                Data = Data/256/numAvg*InputRange + Offset;
            case 5
                Data = loadBinary(Dpath,'int8',0);
                Data = Data/256*InputRange + Offset;
                Data = Data(1:floor(size(Data,1)/256)*256,:);
            case 6
                Data = loadBinary(Dpath,'double',0);
            case 9
                Data = loadBinary(Dpath,'int16',0);
                Data = Data/256/numAvg*InputRange + Offset;
            case 26
                Data = loadBinary(Dpath,'double',0);
            case {31,32}
                Data = loadBinary(Dpath,'double',0);
            case 33
                Data = loadBinary(Dpath,'int16',0);
                Data = Data/256/numAvg*InputRange + Offset;
            case 34
                Data = loadBinary(Dpath,'int8',0);
                Data = Data/256*InputRange + Offset;
                Data = Data(1:floor(size(Data,1)/256)*256,:);
            case 35
                Data = loadBinary(Dpath,'int16',0);
                Data = Data/256/numAvg*InputRange + Offset;
            case 36
                Data = loadBinary(Dpath,'int16',0);
                Data = Data/256/numAvg*InputRange + Offset;
            case 101
                Data = -loadBinary(Dpath,'double',0);
            case default
                display('DataLoader: unknown system ID! Attempting to load data...');
                Data = loadBinary(Dpath,'int16',0);
        end
    end
end
Data = Data'; % transpose data sp the 1st dimension is time
end
