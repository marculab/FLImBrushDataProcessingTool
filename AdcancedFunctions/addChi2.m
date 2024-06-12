clear
close all
clc
addpath(genpath('../Algorithms'))
%% 
root = "G:\V2_Sacramento_Database Apr-30-2024\Da Vinci Robot Study (100 patients)\Data_100_Patient_Study";
patientID = 101:200;
for i = 2:numel(patientID)
    ID = patientID(i);
    folder = dir(fullfile(root,['*' num2str(ID)]));
    DeConFile = dir(fullfile(folder.folder,folder.name,'Triplex/ALL_DECONVOLVED_FILES','*.mat'));
    liteFLag = zeros(size(DeConFile));
    for m = 1:numel(DeConFile)
        nameTemp = DeConFile(m).name;
        liteFLag(m) = contains(nameTemp,'lite');
    end
    DeConFile = DeConFile(~(liteFLag));

    for m = 1:numel(DeConFile)
    input_mat_name = fullfile(DeConFile(m).folder,DeConFile(m).name)
    load(input_mat_name)
    saveDeconLite(Ch1DataObj,Ch2DataObj,Ch3DataObj,'',input_mat_name);
    end
end