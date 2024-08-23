clear
close all
clc

addpath(genpath('..\Algorithms'))
%% load in new APD file
APDFolder = '..\APDDetectorFile\UCD_MC';
apd1Obj = apdClass(fullfile(APDFolder,'M00562768_laser_reflection.mat'));
apd1Obj.creatFromPath();

apd2Obj = apdClass(fullfile(APDFolder,'M00549707_DCS.mat'));
apd2Obj.creatFromPath();

apd3Obj = apdClass(fullfile(APDFolder,'M00549708_DASPI.mat'));
apd3Obj.creatFromPath();

apd4Obj = apdClass(fullfile(APDFolder,'M00928868 DASPI.mat'));
apd4Obj.creatFromPath();
%% load in BG review table
RunLevelData = importHNRunData("H&N Aggregate Input Data 20240822.xlsx", "Run_Level", [1, Inf]);
T = RunLevelData(RunLevelData.Patient>100,:);
%%
root = "F:\V2_Sacramento_Database\Da Vinci Robot Study (100 patients)\Data_100_Patient_Study";
for i = 568%:numel(T.Patient)
    ID = T.Patient(i);
    Run = T.Run(i);
    DeConFile = T.DeconvolutionFile(i)
    DeConFile = fullfile(root,DeConFile);
    ReprocessFB20240822EstimateBG(DeConFile,'A:',apd1Obj,apd2Obj,apd3Obj,apd4Obj) % change the drive letter
    
end