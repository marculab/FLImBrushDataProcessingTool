% code to interplocation of V6 data
%%
close all
clear 
clc

addpath(genpath(pwd))
addpath(genpath('C:\Users\Xiangnan\Documents\MyGitRepo\TRFS_TRIPLEX_GUI\Algorithms'))
%% load in data
root = 'W:\V2 Sacramento Database\Brain Necrosis\2019-FirstStudyDrBloch\RawData\Subject011-20200226'; % defaut folder for file selection
%% set save path
savePath = 'W:\V2 Sacramento Database\Brain Necrosis\2019-FirstStudyDrBloch\RawData\Subject011-20200226\MC Videos';
[DeConfile,DeConpath] = uigetfile([root '\*.mat'],'Please select DeCon file','MultiSelect','off');

fileTemp = DeConfile;
[~,runName,~] = fileparts(fileTemp);
runNum = strsplit(runName,'_');
runNum = str2double(runNum{5});
%% get vidoe file
[vidFile,vidPath] = uigetfile([root '\*.avi'],'Please select video file','MultiSelect','off');

%% get location file
[locFile,locPath] = uigetfile([root '\*.mat'],'Please select location file','MultiSelect','off');
%% get txt file
[txtFile,txtPath] = uigetfile([root '\*.txt'],'Please select txt file','MultiSelect','off');

%% load in data
load(fullfile(DeConpath,fileTemp))
load(fullfile(locPath,locFile))
VideoData = importVideoTextFile(fullfile(txtPath,txtFile));
RepRate = 100;

%% open video and check whether num of frames matches txt file
v = VideoReader(fullfile(vidPath, vidFile));
if v.NumFrames ~= size(VideoData,1)
    errordlg('Number of frames mismatch between video and txt file!','Error');
    error('Number of frames mismatch between video and txt file!');
end
im = read(v,1);
%% match data point to frames
frameIdxVid = 1:size(VideoData,1);
frameIdxVid = frameIdxVid';
VideoData(:,9) = VideoData(:,9)-VideoData(1,9);
frameT = VideoData(:,9);
[frameT,ia] = unique(frameT); % find duplicate frame
frameIdxVid = frameIdxVid(ia); % remove duplicated frame
time = timeStamp*1000; % time vector of data points
FrameIdxData = interp1(frameT,frameIdxVid,time);
FrameIdxData = round(FrameIdxData);

%% loop through all frames and creat overlay image
dataToPlot = lifet_avg{2};
% scale = [floor(quantile(dataToPlot,0.1)) ceil(quantile(dataToPlot,0.9))];
scale = [6 9];
radius = 12;
alpha = 0.5;
numOfFrames = v.NumFrames;

overlayAll = cell(numOfFrames,1);

parfor i = 1:numOfFrames
% use motion correction location data
% dataToPlotTemp = dataToPlot(FrameIdxData<=i);
pos = posMC{i};
if i>3 % if there are more than 3 frames
%------------replace single 0 lines with arerage of before and after---------------------------------------%
for j = 2:i-1
    posCurrent = pos(j,1)+ pos(j,2);
    if posCurrent==0 % if current pos is all 0 check before and after
        posBefore = pos(j-1,1)+ pos(j-1,2);
        posAfter = pos(j+1,1)+ pos(j+1,2);
        if posBefore&&posAfter % if single 0, replace with average
            pos(j,1) = 0.5*(pos(j-1,1)+pos(j+1,1));
            pos(j,2) = 0.5*(pos(j-1,2)+pos(j+1,2));
        end
    end
end
end
% set all 0 points (no aiming beam location) to NaN;
temp = sum(pos,2);
ZeroIdx = find(temp==0);
pos(ZeroIdx,:)=NaN(length(ZeroIdx),2);

%% interpolate locations
% xx = interp1(frameT,pos(:,1),time);
xx = pos(:,1);
% xx(isnan(xx))=0; % replace NaN with 0
% yy = interp1(frameT,pos(:,2),time);
yy = pos(:,2);

% yy(isnan(yy))=0; % replace NaN with 0
% rr(isnan(rr)) = 0;

% set position daya
posData = struct();
posData.px = xx;
posData.py = yy;
posData.frameIdx = FrameIdxData;

% generate overlay 
[overlay,~] = getOverlay(size(im), posData, dataToPlot, scale, radius);
% imshow(overlay)
overlayAll{i} = overlay;
end
%% loop through all frames
cd(savePath)
[filepath,name,ext] = fileparts(fileTemp);
save([name '_MC_overlay_all.mat'],'overlayAll');
%% make video
frameRate = 30;
w = VideoWriter([name '_MC_' num2str(frameRate) 'Hz.mp4'],'MPEG-4');
w.FrameRate = frameRate;
open(w)
%---------------------------------------plot channel 1 video-------------------------------------------------------------------
figure('units','normalized','outerposition',[0.1 0.1 0.8 0.8])

ROI_x = 1:size(im,1);
ROI_y = 1:size(im,2);

plotStep =1;

for i = 1:plotStep:numOfFrames
    im = read(v,i);
    overlayTemp = overlayAll{i};
%     imshow(overlayTemp)
    BW = sum(overlayTemp,3);
    mask= cat(3,BW,BW,BW); 
%     imshow(mask)
    augmentedImgCh1LT = im;
    augmentedImgCh1LT(~(mask == 0)) = alpha*overlayTemp(~(mask == 0) ) + (1-alpha)*im( ~(mask == 0) );
    
    
    % scale = [1 7];
%     [augmentedImgCh1LT,~] = AugmentImg(im, posData, plotTemp, scale, radius, alpha);
    augmentedImgCh1LT = augmentedImgCh1LT(ROI_x,ROI_y,:);
    % % show image
    imshow(augmentedImgCh1LT)
%     title([name ' Channel 1 lifetime'],'Interpreter','none')
    colormap(jet);
    caxis(gca,scale);
    h0 = colorbar;
    % %     ylabel(h0, ['Lifetime CH', int2str(dest_channel),' (ns)'])
    h0.Label.String = 'Lifetime (ns)';
    set(gca,'FontSize',12)
    set(gca,'LooseInset',get(gca,'TightInset'))
    frame = getframe(gcf);
    writeVideo(w,frame);
end
close(w)

%%
cd(root)
