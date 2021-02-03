% get white light image from video
clear all
close all
clc

%% get filels
[file,path] = uigetfile('.avi','MultiSelect','on');

%% open video and save image
cd(path)
for i = 1:length(file)
    name = file{i};
    fullname = fullfile(path,name);
    [pathstr, saveName, ext] = fileparts(name);
    v = VideoReader(fullname);
    D= v.Duration-0.040;
    v.CurrentTime = D;
%     v = VideoReader(name,'CurrentTime',D);
    im = readFrame(v);
    imshow(im)
    drawnow
    print([saveName '.jpg'],'-djpeg','-r600')
end