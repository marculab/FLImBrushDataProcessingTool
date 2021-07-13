%%
close all
clear all
clc

addpath(genpath(pwd))

%% load in data
root = 'C:\Users\Xiangnan\Box Sync\V6 Clinical Data\20210505_neuro024';
DeConpath = uigetdir(root,'Please select reconstructed .mat file');
VidPath = uigetdir(root,'Please select video file');
files = dir(fullfile(DeConpath,'*ImgRecon.mat'));

vidFiles = dir(fullfile(VidPath,'*.avi'));
TF = false(size(vidFiles));
for i = 1:size(vidFiles,1)
    TF(i) = ~contains(vidFiles(i).name,'augmented');
end
vidFiles = vidFiles(TF);
%%
cd(DeConpath)
% replotVideo(['videos\' videoName '.avi'], [name '_interp.mat'])
for i = 1:length(files)
    tempName = files(i).name;
    scale = cell(3,1);
    scale{1} = [2 5];
    scale{2} = [2 5];
    scale{3} = [2 5];
    [filepath,name,ext] = fileparts(tempName);
    temp = split(name,'_');
    VidName = vidFiles(i).name;
    replotInterpImage([name '.mat'],fullfile(VidPath,VidName),10,0.5,20,scale)
end
close all