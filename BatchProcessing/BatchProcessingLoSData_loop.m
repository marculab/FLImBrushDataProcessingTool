clear
close all
clc
addpath(genpath('..\Algorithms'))

APD1Path = '..\APDDetectorFile\UCD_MC\M00562768_laser_reflection.mat';
APD2Path = '..\APDDetectorFile\UCD_MC\M00549707_DCS.mat';
APD3Path = '..\APDDetectorFile\UCD_MC\M00549708_DASPI.mat';
APD4Path = '..\APDDetectorFile\UCD_MC\M00928868 DASPI.mat';

%Decon variables
ChWidth = 680*0.08;
alpha1 = 0.916;
alpha2 = 0.916;
alpha3 = 0.916;
alpha4 = 0.916;
%LGorder = 12;

%loop variables 
% TBin = [680, 800, 925, 1050, 1175, 1300, 1425, 1550, 1675, 1800, 1925];
% ChWidth = TBin *0.08;
% alpha = [0.80,0.82,0.84,0.86,0.88,0.90,0.916,0.93,0.95,0.965,0.98,0.99];
LGorder = [6, 8, 10, 12, 14, 16, 18];

root = 'C:\Users\MarcuAdmin\Documents\Data\LoS_2024\240710-C440+C460\FLImBrush\';
Datafolder = 'C:\Users\MarcuAdmin\Documents\Data\LoS_2024\240710-C440+C460\FLImBrush';
BGFile = 'C:\Users\MarcuAdmin\Documents\Data\LoS_2024\240710-C440+C460\FLImBrush\Background\Background01.tdms';

dataFiles = dir([Datafolder '\*.tdms']);
%% 
%for i = 4:numel(dataFiles)
for i = 16:16

    disp(['Processing ' dataFiles(i).name]);
    DataPath = fullfile(dataFiles(i).folder, dataFiles(i).name);

    % for j = 1:numel(ChWidth)
    % ChWidth_j = ChWidth(j);
    % processRun_HN(root, DataPath, BGFile, APD1Path, APD2Path, APD3Path, APD4Path, alpha1, alpha2, alpha3, alpha4, ChWidth_j, LGorder);
    % end

    % for j = 1:numel(alpha)
    % alpha_j = alpha(j)
    % processRun_HN(root, DataPath, BGFile, APD1Path, APD2Path, APD3Path, APD4Path, alpha_j, alpha_j, alpha_j, alpha_j, ChWidth, LGorder);
    % end

    for j = 1:numel(LGorder)
    LGorder_j = LGorder(j);
    processRun_HN(root, DataPath, BGFile, APD1Path, APD2Path, APD3Path, APD4Path, alpha1, alpha2, alpha3, alpha4, ChWidth, LGorder_j);
    end

end