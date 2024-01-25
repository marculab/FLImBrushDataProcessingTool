clear
close all
clc

saveRootPath = 'C:\Users\MarcuAdmin\Documents\Data\BatchProcess';
DataFilePath = 'F:\V2_Sacramento_Database\Da Vinci Robot Study (100 patients)\Data_100_Patient_Study\';
BGFilePath = 'F:\V2_Sacramento_Database\Da Vinci Robot Study (100 patients)\Data_100_Patient_Study\';
APD1Path = 'C:\Users\MarcuAdmin\Documents\GitHub\FLImBrushDataProcessingTool\APDDetectorFile\UCD_MC\M00562768_laser_reflection.mat';
APD2Path = 'C:\Users\MarcuAdmin\Documents\GitHub\FLImBrushDataProcessingTool\APDDetectorFile\UCD_MC\M00549707_DCS.mat';
APD3Path = 'C:\Users\MarcuAdmin\Documents\GitHub\FLImBrushDataProcessingTool\APDDetectorFile\UCD_MC\M00549708_DASPI.mat';
ChWidth = 54.4;
alpha1 = 0.916;
alpha2 = 0.916;
alpha3 = 0.916;

DataAggregate = importDataAggregate('H&N Aggregate Input Data - Run_Level.csv');
% DataAggregate(1,:)=[]; 
%% 

serverDriveLetter = 'C:';
%for i = 1:size(DataAggregate,1)
for i = 1009:1068
%    if DataAggregate.RawFLImDataFiletdms{i} ~= "N/A"
    if ~isempty(DataAggregate.RawFLImDataFiletdms{i}) && ~strcmpi(DataAggregate.RawFLImDataFiletdms{i}, "N/A")
        DataPath = [DataFilePath, DataAggregate.RawFLImDataFiletdms{i}];
        disp(['Processing ' DataPath]);
        BGPath = [BGFilePath, char(DataAggregate.BackgroundFile(i))];
        % ChannelUsed = string(DataAggregate.DataChannelsUsed(i));
        % if contains(ChannelUsed,'CH4')
        %     alpha3 = 0.916;
        % else
        %     alpha3 = 0.916;
        % end
        processRun_HN_KE(DataPath, BGPath, APD1Path, APD2Path, APD3Path, alpha1, alpha2, alpha3, ChWidth, saveRootPath);
        clear DataPath BGPath;

    end
end