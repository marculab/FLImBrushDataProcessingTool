clear
close all
clc
addpath(genpath('..\Algorithms'))
saveRootPath = 'C:\Users\Xiangnan\Box\Prostate Pilot\ProstateBxDataBase';
% DataPath = 'Z:\V2 Sacramento Database\SP4-Intuitive Surgical\2023-03-21-Visit\2023_03_22_PEM\2023_03_22_PEM_50.tdms';
% BGPath = 'Z:\V2 Sacramento Database\SP4-Intuitive Surgical\2023-03-21-Visit\2023_03_21_Baseline_Tissue\Background\Background01.tdms';
APD1Path = '..\APDDetectorFile\UCD_NCIBT\M00936636.mat';
APD2Path = '..\APDDetectorFile\UCD_NCIBT\M00972128.mat';
APD3Path = '..\APDDetectorFile\UCD_NCIBT\M00975073.mat';
APD4Path = '..\APDDetectorFile\UCD_NCIBT\M00936637.mat';
ChWidth = 680*0.08;
alpha1 = 0.916;
alpha2 = 0.916;
alpha3 = 0.916;
alpha4 = 0.916;
Datafolder = 'C:\Users\Xiangnan\Box\Prostate Pilot\ProstateBxDataBase\20240604P4';
BGFile = 'C:\Users\Xiangnan\Box\Prostate Pilot\ProstateBxDataBase\20240604P1\Background\Background01.tdms';
% IntuitiveDataAggregate(1,:)=[];
dataFiles = dir([Datafolder '\*.tdms']);
for i = 1:numel(dataFiles)

    disp(['Processing ' dataFiles(i).name]);
    DataPath = fullfile(dataFiles(i).folder, dataFiles(i).name);
    processRun(DataPath, BGFile, APD1Path, APD2Path, APD3Path, APD4Path, alpha1, alpha2, alpha3, alpha4, ChWidth, saveRootPath);

end