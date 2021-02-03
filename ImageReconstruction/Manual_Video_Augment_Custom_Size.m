clear all;
close all;
home;
addpath(genpath(pwd))
%% IMPORT VIDEO FILE
hold on
obj=VideoReader('20200908LampChop0.24mW_run20.avi'); % Specify the video file to load
I=read(obj); % Stores data in variable I
% implay(I,18); 
video_size=size(I);
max_frame=video_size(4);
for n=1:max_frame
    figure(1);
    frame_temp=I(:,:,:,n);
    frame_temp=imresize(frame_temp,[720 1280]); % Enter Height & Width of Video Dimensions e.g. 720 x 1280
    imshow(frame_temp);
    [x,y] = ginput(1)
    A=x>0;
    B=x<1280;
    C=y>0;
    D=y<720;
    Filt=A.*B.*C.*D;
    x=Filt.*x;
    y=Filt.*y;
    X(n)=x;
    Y(n)=y;
end
Coordinates=[X; Y]' % This is the output file of interest
%% Extract Video Frame
obj=VideoReader('P528731_S_07_27_2020_run7.avi'); % Specify the video file to load
I=read(obj); % Stores data in variable I
implay(I,18); 
%% Write Video Frame
thisFrame = read(obj,7847); % Enter Frame Number on Right
thisFrame=imresize(thisFrame,[720 1280]); % Resize Frame if not 720x1080
imwrite(thisFrame,'Run_7_Raw.jpeg')