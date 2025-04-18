clear
close all
clc

addpath(genpath('E:\MyGitRepo\FLImBrushDataProcessingTool\Algorithms'))
%% load in new APD file
APDFolder = 'C:\Users\Xiangnan\Box\APD Detector data file backup\UCD MC\20240322 APD Files after 4Ch upgrade';
apd1Obj = apdClass(fullfile(APDFolder,'M00562768_laser_reflection 03-21-2024 14-11.mat'));
apd1Obj.creatFromPath();

apd2Obj = apdClass(fullfile(APDFolder,'M00549707_DCS 03-21-2024 14-53.mat'));
apd2Obj.creatFromPath();

apd3Obj = apdClass(fullfile(APDFolder,'M00549708_DASPI 03-21-2024 15-15.mat'));
apd3Obj.creatFromPath();
%% load in BG review table
BkgdreviewApril16 = importBGReview("Bkgd-review-April_16.xlsx");
T = BkgdreviewApril16(BkgdreviewApril16.selectNewBkgd=="x",:);
%%
root = "F:\V2_Sacramento_Database\Da Vinci Robot Study (100 patients)\Data_100_Patient_Study";
for i = 1:numel(T.Patient)
    ID = T.Patient(i);
    Run = T.Run(i);
    folder = dir(fullfile(root,['*' num2str(ID)]));
    DeConFile = dir(fullfile(folder.folder,folder.name,'Triplex/ALL_DECONVOLVED_FILES',sprintf('*_%02d_12.5GS.mat',Run)));
    for m = 1:numel(DeConFile)
    input_mat_name = fullfile(DeConFile(m).folder,DeConFile(m).name)
    ReprocessFB20240426EstimateBG(input_mat_name,apd1Obj,apd2Obj,apd3Obj)
    end
end