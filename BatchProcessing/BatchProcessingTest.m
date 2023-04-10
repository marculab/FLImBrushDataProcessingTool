clear
close all
clc

saveRootPath = 'C:\Users\Xiangnan\Desktop\BatchProcessIntuitive';
% DataPath = 'Z:\V2 Sacramento Database\SP4-Intuitive Surgical\2023-03-21-Visit\2023_03_22_PEM\2023_03_22_PEM_50.tdms';
% BGPath = 'Z:\V2 Sacramento Database\SP4-Intuitive Surgical\2023-03-21-Visit\2023_03_21_Baseline_Tissue\Background\Background01.tdms';
APD1Path = 'E:\MyGitRepo\FLImBrushDataProcessingTool\APDDetectorFile\UCD_GBSF\M00570374.mat';
APD2Path = 'E:\MyGitRepo\FLImBrushDataProcessingTool\APDDetectorFile\UCD_GBSF\M00570375.mat';
APD3Path = 'E:\MyGitRepo\FLImBrushDataProcessingTool\APDDetectorFile\UCD_GBSF\M00570376.mat';
ChWidth = 154;
alpha1 = 0.916;
alpha2 = 0.916;
% alpha3 = 0.916;

IntuitiveDataAggregate = importIntuitiveDataAggregate('Intuitive Surgical Data Aggregate.csv');
% IntuitiveDataAggregate(1,:)=[];
serverDriveLetter = 'Z:';
for i = 84:size(IntuitiveDataAggregate,1)
    if IntuitiveDataAggregate.RawDataFile{i} ~= "N/A"
        DataPath = [serverDriveLetter, IntuitiveDataAggregate.RawDataFile{i}];
        disp(['Processing ' DataPath]);
        BGPath = [serverDriveLetter, char(IntuitiveDataAggregate.BackgroundFile(i))];
        ChannelUsed = string(IntuitiveDataAggregate.DataChannelsUsed(i));
        if contains(ChannelUsed,'CH4')
            alpha3 = 0.965;
        else
            alpha3 = 0.916;
        end
        processRun(DataPath, BGPath, APD1Path, APD2Path, APD3Path, alpha1, alpha2, alpha3, ChWidth, saveRootPath);
    end
end