clear
close all
clc

addpath(genpath('..\Algorithms'))
%% load in new APD file
APDFolder = '..\APDDetectorFile\UCD_MC';
apd1Obj = apdClass(fullfile(APDFolder,'M00562768_LASER.mat'));
apd1Obj.creatFromPath();

apd2Obj = apdClass(fullfile(APDFolder,'M00549707_DCS.mat'));
apd2Obj.creatFromPath();

apd3Obj = apdClass(fullfile(APDFolder,'M00549708_DASPI.mat'));
apd3Obj.creatFromPath();

apd4Obj = apdClass(fullfile(APDFolder,'M00928868 DASPI.mat'));
apd4Obj.creatFromPath();
%% load in BG review table
RunLevelData = importHNRunData("H&N Aggregate Input Data 20241029.xlsx", "Run_Level", [1, Inf]);
T = RunLevelData(RunLevelData.Patient>100,:);
%%
root = "A:\V2_Sacramento_Database\Da Vinci Robot Study (100 patients)\Data_100_Patient_Study";
for i = 575:575
    ID = T.Patient(i);
    Run = T.Run(i);
    DeConFile = T.DeconvolutionFile(i);
    if ~(DeConFile=="")
        if isfile(fullfile(root,DeConFile))
            DeConFile = fullfile(root,DeConFile);
        else
            [filepath,name,ext] = fileparts(DeConFile);
            temp = dir(fullfile(root,filepath, strcat(name, '*')))
            Idx = 0;
            for m = 1:length(temp)
                if ~contains(temp(m).name,'lite')
                    Idx = m;
                end
            end
            DeConFile = fullfile(root,filepath,temp(Idx).name);
        end
        DeConFile
        expOrder = 3; % exp decon order
        tic
        ReprocessFB20241030EstimateBG(DeConFile,'A:',apd1Obj,apd2Obj,apd3Obj,apd4Obj,expOrder) % change the drive letter
        toc
    end
end