clear all;
close all;
home;
addpath(genpath(pwd))
%% IMPORT VIDEO FILE
[file,path] = uigetfile('*.avi','Please select avi video file');
obj=VideoReader(fullfile(path,file)); % Specify the video file to load
%% read in video data
I=read(obj); % Stores data in variable I
display('Finished reading video file!')
video_size=size(I);
display(video_size);
max_frame=video_size(4);
%% load in CNN based segmentation location
[file,path] = uigetfile('*.mat','Please select CNN segmentation file');
load(fullfile(path,file))
pos = pos(:,7:8);
pos = double(pos);
display('Finished loading CNN result!')

%% find frame with bad segmentation
ZeroIdx = find(~sum(pos,2));
%% show CNN segmentation result
fig = figure('Position',[200 320 1280 720]);

for n=1:max_frame
    frame_temp=I(:,:,:,n);
    frame_temp=imresize(frame_temp,[720 1280]); % Enter Height & Width of Video Dimensions e.g. 720 x 1280
    imshow(frame_temp,'Parent',gca);
    viscircles(pos(n,:),2)
    title(['Frame ' num2str(n)])
    drawnow limitrate
    
    %     [x,y] = ginput(1)
    %     A=x>0;
    %     B=x<1280;
    %     C=y>0;
    %     D=y<720;
    %     Filt=A.*B.*C.*D;
    %     x=Filt.*x;
    %     y=Filt.*y;
    %     X(n)=x;
    %     Y(n)=y;
end
% Coordinates=[X; Y]' % This is the output file of interest
%% correct bad frames
X = zeros(size(pos,1),1);
Y = X;
for n= 1:size(pos,1)
    frameNum = n;
    frame_temp=I(:,:,:,frameNum);
    frame_temp=imresize(frame_temp,[720 1280]); % Enter Height & Width of Video Dimensions e.g. 720 x 1280
    imshow(frame_temp,'Parent',gca);
    viscircles(pos(n,:),2)
    title(['Frame ' num2str(frameNum)])
    drawnow limitrate
    [x,y] = ginput(1);
    A=x>0;
    B=x<1280;
    C=y>0;
    D=y<720;
    Filt=A.*B.*C.*D;
    if Filt
        X(n)=x;
        Y(n)=y;
    else
        X(n)=pos(n,1);
        Y(n)=pos(n,2);
    end
end
pos_new=[X Y];

