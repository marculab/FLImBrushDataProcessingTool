% code to interplocation of V6 data
%%
close all
clear 
clc

addpath(genpath(pwd))
addpath(genpath('..'))
%% load in data
root = 'D:\LoaclData\5ALAFLImTest\Subject031_20210804'; % defaut folder for file selection
%% set save path
savePath = 'C:\Users\Xiangnan\Desktop\OverlayTest';

[DeConfile,DeConpath] = uigetfile([root '\*.mat'],'Please select DeCon file','MultiSelect','off');


fileTemp = DeConfile;
[~,runName,~] = fileparts(fileTemp);
runNum = strsplit(runName,'_');
runNum = str2double(runNum{5});

plotCh1 = 1;
plotCh2 = 1;
plotCh3 = 1;
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
RepRate = Ch1DataObj.laserRepRate;

%% open video and check whether num of frames matches txt file
v = VideoReader(fullfile(vidPath, vidFile));
if v.NumFrames ~= size(VideoData,1)
    errordlg('Number of frames mismatch between video and txt file!','Error');
    error('Number of frames mismatch between video and txt file!');
end
im = read(v,1);
%% load in fluorescence PPIX image
[fmFile,fmPath] = uigetfile([root '\*'],'Please select fluorescence microscope image.','MultiSelect','off');
%%
FMRaw = imread(fullfile(fmPath,fmFile));
FMRaw = FMRaw(445:1092,:,:);
figure;imshow(FMRaw);
F = griddedInterpolant(double(FMRaw));
[sx,sy,sz] = size(FMRaw);
xq = (0:1.35:sx-1.35)';
yq = (0:1.35:sy-1.35)';
zq = (1:sz)';
FM = uint8(F({xq,yq,zq}));
%% calculate num of data points
% shift = length(Ch1LT)-(MetaData(end,9)-MetaData(1,9))/1000*120/4;
% shift = round(shift);
shift=0; %shift need to be 0!

Ch1INTCorr = circshift(Ch1DataObj.INTsAllGainCorrected,shift);
Ch1LT = circshift(Ch1DataObj.LTsAll,shift);
Ch1SNR = circshift(Ch1DataObj.SNR,shift)';
G1 = circshift(Ch1DataObj.gain,shift);

Ch2INTCorr = circshift(Ch2DataObj.INTsAllGainCorrected,shift);
Ch2LT = circshift(Ch2DataObj.LTsAll,shift);
Ch2SNR = circshift(Ch2DataObj.SNR,shift)';
G2 = circshift(Ch2DataObj.gain,shift);

Ch3INTCorr = circshift(Ch3DataObj.INTsAllGainCorrected,shift);
Ch3LT = circshift(Ch3DataObj.LTsAll,shift);
Ch3SNR = circshift(Ch3DataObj.SNR,shift)';
G3 = circshift(Ch3DataObj.gain,shift);

%% filter data by gain
G1Max = 250; % due to BG subtraction
G2Max = 2000;
G3Max = 2000;

Ch1INTCorr(G1>G1Max) = NaN;
Ch1LT(G1>G1Max) = NaN;
Ch1SNR(G1>G1Max) = NaN;

Ch2INTCorr(G2>G2Max) = NaN;
Ch2LT(G2>G2Max) = NaN;
Ch2SNR(G2>G2Max) = NaN;

Ch3INTCorr(G3>G3Max) = NaN;
Ch3LT(G3>G3Max) = NaN;
Ch3SNR(G3>G3Max) = NaN;
%% match data point to frames
frameIdxVid = 1:size(VideoData,1);
frameIdxVid = frameIdxVid';
VideoData(:,9) = VideoData(:,9)-VideoData(1,9);
frameT = VideoData(:,9);
[frameT,ia] = unique(frameT); % find duplicate frame
frameIdxVid = frameIdxVid(ia); % remove duplicated frame
time = 0:1:size(Ch1LT,1)-1; % time vector of data points
time = time'*1000/RepRate*4; % time vector of data points
FrameIdxData = interp1(frameT,frameIdxVid,time);
FrameIdxData = round(FrameIdxData);

%% loop through all frames and creat overlay image
dataToPlot = Ch3LT;
% scale = [floor(quantile(dataToPlot,0.1)) ceil(quantile(dataToPlot,0.9))];
scale = [0 12];
radius = 7.5;
alpha = 0.5;
numOfFrames = v.NumFrames;

overlayAll = cell(numOfFrames,1);

plotStep =1;
parfor i = 1:1:numOfFrames
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
xx = interp1(frameT,pos(:,1),time);
% xx(isnan(xx))=0; % replace NaN with 0
yy = interp1(frameT,pos(:,2),time);
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
frameRate = 60;
w = VideoWriter([name '_MC_' num2str(frameRate) 'Hz.mp4'],'MPEG-4');
w.FrameRate = frameRate;
open(w)
%---------------------------------------plot channel 1 video-------------------------------------------------------------------
figure('units','normalized','outerposition',[0.1 0.1 0.8 0.8])

ROI_x = 1:size(im,1);
ROI_y = 1:size(im,2);

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
overlayTemp = overlayAll{1701};
BW = sum(overlayTemp,3);
mask= cat(3,BW,BW,BW);
augmentedFMImgCh1LT = FM;
augmentedFMImgCh1LT(~(mask == 0)) = alpha*overlayTemp(~(mask == 0) ) + (1-alpha)*augmentedFMImgCh1LT( ~(mask == 0) );
augmentedFMImgCh1LT = augmentedFMImgCh1LT(ROI_x,ROI_y,:);
% % show image
imshow(augmentedFMImgCh1LT)
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
close(w)

%%
cd(root)
