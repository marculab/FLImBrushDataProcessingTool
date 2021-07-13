clear all;
close all;
home;
addpath(genpath(pwd))
%% IMPORT VIDEO FILE
root = 'C:\Users\Xiangnan\Box Sync\FLImBrush vs V4\20210505';
[file,path] = uigetfile([root '\*.avi'],'Please select avi video file');
obj=VideoReader(fullfile(path,file)); % Specify the video file to load
%% read in video data
I=read(obj); % Stores data in variable I
display('Finished reading video file!')
video_size=size(I);
display(video_size);
max_frame=video_size(4);
%% select ROI
figure
frame_temp=I(:,:,:,1);
frame_temp=imresize(frame_temp,[obj.Height obj.Width]); % Enter Height & Width of Video Dimensions e.g. 720 x 1280
imshow(frame_temp,'Parent',gca);
ROI = drawfreehand(gca);
%% load in CNN based segmentation location
[file,path] = uigetfile([root '\*.mat'],'Please select CNN segmentation file');
load(fullfile(path,file))
pos = pos(:,7:8);
pos = double(pos);
display('Finished loading CNN result!')

%% find frame with bad segmentation
ZeroIdx = find(~sum(pos,2));
%% show CNN segmentation result
fig1 = figure('Position',[200 320 obj.Height obj.Width]);

for n=1:max_frame
    frame_temp=I(:,:,:,n);
    frame_temp=imresize(frame_temp,[obj.Height obj.Width]); % Enter Height & Width of Video Dimensions e.g. 720 x 1280
    imshow(frame_temp,'Parent',gca);
    viscircles(pos(n,:),2)
    title(['Frame ' num2str(n) ' x ' num2str(pos(n,1)) ' y ' num2str(pos(n,2))])
    drawnow limitrate
end
close(fig1)
% Coordinates=[X; Y]' % This is the output file of interest
%% correct bad frames
X = zeros(size(pos,1),1);
Y = X;
fig2 = figure('Position',[200 320 obj.Height obj.Width]);
for n= 1:size(pos,1)
    frameNum = n;
    frame_temp=I(:,:,:,frameNum);
    frame_temp=imresize(frame_temp,[obj.Height obj.Width]); % Enter Height & Width of Video Dimensions e.g. 720 x 1280
    imshow(frame_temp,'Parent',gca);
    viscircles(pos(n,:),2)
    title(['Frame ' num2str(n) ' x ' num2str(pos(n,1)) ' y ' num2str(pos(n,2))])
    drawnow limitrate
    [x,y] = ginput(1);
    A=x>0;
    B=x<obj.Width;
    C=y>0;
    D=y<obj.Height;
    Filt=A.*B.*C.*D;
    if Filt
        X(n)=x;
        Y(n)=y;
    else
        X(n)=pos(n,1);
        Y(n)=pos(n,2);
    end
end
close(fig2)
pos_new=[X Y];


%% remove points out of ROI
TF = inROI(ROI,pos_new(:,1),pos_new(:,2));
TF = ~TF;
pos_new(TF,:) = zeros(sum(TF),2);
save(fullfile(path,'pos_corrected'),'pos_new')



