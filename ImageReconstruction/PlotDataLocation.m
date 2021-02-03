% code to plot data location as white dot
%%
clear all
close all
clc
addpath(genpath(pwd))

%% selection data and video file
root = 'C:\Users\Xiangnan\OneDrive - University of California, Davis\UC Davis\Paper&Grant\My Paper\2020 APD Paper\Data';
[dataFile,dataPath] = uigetfile([root '\*.mat'],'Please select DeCon file');
[vidFile,vidPath] = uigetfile([root '\*.avi'],'Please select video file');

%% load image
v = VideoReader(fullfile(vidPath,vidFile));
D= v.Duration-0.04;
v.CurrentTime = D;
im = readFrame(v);
imshow(im)
drawnow
%% load location
load(fullfile(dataPath,dataFile))
% outputMat = outputMat(1:20:end,:);
posData = {};
posData.px = outputMat(:, 7);
posData.py = outputMat(:, 8);
posData.frames = outputMat(:,1);
posData.radius = outputMat(:,9)*0.1;

ltData = {};
ltData.lt  = {outputMat(:, 12), outputMat(:, 13), outputMat(:, 14), outputMat(:, 15)};
ltData.it  = {outputMat(:, 16), outputMat(:, 17), outputMat(:, 18), outputMat(:, 19)};
ltData.snr = {outputMat(:, 20), outputMat(:, 21), outputMat(:, 22), outputMat(:, 23)};
cd(vidPath)
[filepath,name,ext] = fileparts(dataFile);
processImgLocation(im, posData, ltData, name);
