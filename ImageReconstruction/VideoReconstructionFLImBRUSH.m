% code to interplocation of V6 data
%%
close all
clear 
clc

addpath(genpath(pwd))
% addpath(genpagenpth('..\'))
%% load in data
root = 'D:\LoaclData\5ALAFLImTest\Subject032_20210811'; % defaut folder for file selection
%% set save path
savePath = 'D:\LoaclData\5ALAFLImTest\Subject032_20210811\_OverlayVideoL';

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
figure;imshow(FM);

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

%% repopulate output data mat with interplation
% reset start time to 0
VideoData(:,9) = VideoData(:,9)-VideoData(1,9);
frameTime = VideoData(:,9);
% use CNN segmentation location data
try
    VideoData(:,6:7) = double(pos(:,7:8));
catch
    VideoData(:,6:7) = double(pos_new(:,1:2));
end
% Length = 70000;
% Length = VideoData(end,9);
% VideoData(VideoData(:,9)>Length,:) = [];
numOfVideoFrame = size(VideoData,1);
%------------replace single 0 lines with arerage of before and after---------------------------------------%
for i = 2:numOfVideoFrame-1
    posCurrent = VideoData(i,6)+ VideoData(i,7);
    if posCurrent==0 % if current pos is all 0 check before and after
        posBefore = VideoData(i-1,6)+ VideoData(i-1,7);
        posAfter = VideoData(i+1,6)+ VideoData(i+1,7);
        if posBefore&&posAfter % if single 0, replace with average
            VideoData(i,6) = 0.5*(VideoData(i-1,6)+VideoData(i+1,6));
            VideoData(i,7) = 0.5*(VideoData(i-1,7)+VideoData(i+1,7));
        end
    end
end

% set all 0 points (no aiming beam location) to NaN;
temp = VideoData(:,6)+VideoData(:,7);
temp(temp == min(temp))=0; %set min value to 0 to remove crop original
ZeroIdx = find(temp==0);
VideoData(ZeroIdx,6:7)=NaN(length(ZeroIdx),2);

frameNum = 1:size(VideoData,1);
frameNum = frameNum';

time = 0:1:size(Ch1LT,1)-1;
time = time'*1000/RepRate*4;
frameT = VideoData(:,9);
[frameT,ia] = unique(frameT); % find duplicate frame
VideoData = VideoData(ia,:); % remove duplicated data
frameNum = frameNum(ia); % remove duplicated frame
frameIdx = interp1(frameT,frameNum,time);
frameIdx = round(frameIdx);
output.frame = frameIdx;
VideoData = [frameNum VideoData];


%% interpolate locations
xx = interp1(VideoData(:,10),VideoData(:,7),time);
% xx(isnan(xx))=0; % replace NaN with 0
yy = interp1(VideoData(:,10),VideoData(:,8),time);
% yy(isnan(yy))=0; % replace NaN with 0
rr = interp1(VideoData(:,10),VideoData(:,9),time);
% rr(isnan(rr)) = 0;
output.xx = xx;
output.yy = yy;
output.rr = rr;
output.lt1 = Ch1LT;
output.int1 = Ch1INTCorr;
output.lt2 = Ch2LT;
output.int2 = Ch2INTCorr;
output.lt3 = Ch3LT;
output.int3 = Ch3INTCorr;
output.snr1 = Ch1SNR';
output.snr2 = Ch2SNR';
output.snr3 = Ch3SNR';
output.gain1 = G1;
output.gain2 = G2;
output.gain3 = G3;

%% set position daya
posData.px = xx;
posData.py = yy;
posData.frameIdx = frameIdx;

%% generate overlay for all channels
radius = 7.5;
alpha = 0.5;
scale1 = [floor(quantile(Ch1LT,0.1)) ceil(quantile(Ch1LT,0.9))];
[overlayCh1,overlayCh1All] = getOverlay(size(im), posData, Ch1LT, scale1, radius);

scale2 = [floor(quantile(Ch2LT,0.1)) ceil(quantile(Ch2LT,0.9))];
[overlayCh2,overlayCh2All] = getOverlay(size(im), posData, Ch2LT, scale2, radius);

% scale3 = [floor(quantile(Ch3LT,0.1)) ceil(quantile(Ch3LT,0.9))];
scale3 = [0 20];
[overlayCh3,overlayCh3All] = getOverlay(size(im), posData, Ch3LT, scale3, radius);

ROI_x = 1:size(im,1);
ROI_y = 1:size(im,2);

%% loop through all frames
cd(savePath)
[filepath,name,ext] = fileparts(fileTemp);
%% channel 1
if plotCh1
w = VideoWriter([name '_Ch1.mp4'],'MPEG-4');

open(w)
%---------------------------------------plot channel 1 video-------------------------------------------------------------------
figure('units','normalized','outerposition',[0.1 0.1 0.8 0.8])


for i = 1:1:v.NumFrames
    im = read(v,i);
    overlayTemp = overlayCh1All{i};
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
    caxis(gca,scale1);
    h0 = colorbar;
    % %     ylabel(h0, ['Lifetime CH', int2str(dest_channel),' (ns)'])
    h0.Label.String = 'Lifetime (ns)';
    set(gca,'FontSize',12)
    set(gca,'LooseInset',get(gca,'TightInset'))
    frame = getframe(gcf);
    writeVideo(w,frame);
    
end
overlayTemp = overlayCh1All{end};
BW = sum(overlayTemp,3);
mask= cat(3,BW,BW,BW);
augmentedFMImgCh1LT = FM;
augmentedFMImgCh1LT(~(mask == 0)) = alpha*overlayTemp(~(mask == 0) ) + (1-alpha)*augmentedFMImgCh1LT( ~(mask == 0) );
augmentedFMImgCh1LT = augmentedFMImgCh1LT(ROI_x,ROI_y,:);
% % show image
imshow(augmentedFMImgCh1LT)
%     title([name ' Channel 1 lifetime'],'Interpreter','none')
colormap(jet);
caxis(gca,scale1);
h0 = colorbar;
% %     ylabel(h0, ['Lifetime CH', int2str(dest_channel),' (ns)'])
h0.Label.String = 'Lifetime (ns)';
set(gca,'FontSize',12)
set(gca,'LooseInset',get(gca,'TightInset'))
frame = getframe(gcf);
writeVideo(w,frame);
close(w)
end
%% channel 2
if plotCh2
w = VideoWriter([name '_Ch2.mp4'],'MPEG-4');
open(w)
%---------------------------------------plot channel 1 video-------------------------------------------------------------------
figure('units','normalized','outerposition',[0.1 0.1 0.8 0.8])

for i = 1:1:v.NumFrames
    im = read(v,i);
    overlayTemp = overlayCh2All{i};
%     imshow(overlayTemp)
    BW = sum(overlayTemp,3);
    mask= cat(3,BW,BW,BW); 
%     imshow(mask)
    augmentedImgCh2LT = im;
    augmentedImgCh2LT(~(mask == 0)) = alpha*overlayTemp(~(mask == 0) ) + (1-alpha)*im( ~(mask == 0) );
    
    
    % scale = [1 7];
%     [augmentedImgCh1LT,~] = AugmentImg(im, posData, plotTemp, scale, radius, alpha);
    augmentedImgCh2LT = augmentedImgCh2LT(ROI_x,ROI_y,:);
    % % show image
    imshow(augmentedImgCh2LT)
%     title([name ' Channel 1 lifetime'],'Interpreter','none')
    colormap(jet);
    caxis(gca,scale2);
    h0 = colorbar;
    % %     ylabel(h0, ['Lifetime CH', int2str(dest_channel),' (ns)'])
    h0.Label.String = 'Lifetime (ns)';
    set(gca,'FontSize',12)
    set(gca,'LooseInset',get(gca,'TightInset'))
    frame = getframe(gcf);
    writeVideo(w,frame);
    
end
overlayTemp = overlayCh2All{end};
BW = sum(overlayTemp,3);
mask= cat(3,BW,BW,BW);
augmentedFMImgCh2LT = FM;
augmentedFMImgCh2LT(~(mask == 0)) = alpha*overlayTemp(~(mask == 0) ) + (1-alpha)*augmentedFMImgCh2LT( ~(mask == 0) );
augmentedFMImgCh2LT = augmentedFMImgCh2LT(ROI_x,ROI_y,:);
% % show image
imshow(augmentedFMImgCh2LT)
%     title([name ' Channel 1 lifetime'],'Interpreter','none')
colormap(jet);
caxis(gca,scale2);
h0 = colorbar;
% %     ylabel(h0, ['Lifetime CH', int2str(dest_channel),' (ns)'])
h0.Label.String = 'Lifetime (ns)';
set(gca,'FontSize',12)
set(gca,'LooseInset',get(gca,'TightInset'))
frame = getframe(gcf);
writeVideo(w,frame);
close(w)
end
%% channel 3
if plotCh3
w = VideoWriter([name '_Ch3.mp4'],'MPEG-4');
open(w)
%---------------------------------------plot channel 1 video-------------------------------------------------------------------
figure('units','normalized','outerposition',[0.1 0.1 0.8 0.8])

for i = 1:1:v.NumFrames
    im = read(v,i);
    overlayTemp = overlayCh3All{i};
%     imshow(overlayTemp)
    BW = sum(overlayTemp,3);
    mask= cat(3,BW,BW,BW); 
%     imshow(mask)
    augmentedImgCh3LT = im;
    augmentedImgCh3LT(~(mask == 0)) = alpha*overlayTemp(~(mask == 0) ) + (1-alpha)*im( ~(mask == 0) );
    
    
    % scale = [1 7];
%     [augmentedImgCh1LT,~] = AugmentImg(im, posData, plotTemp, scale, radius, alpha);
    augmentedImgCh3LT = augmentedImgCh3LT(ROI_x,ROI_y,:);
    % % show image
    imshow(augmentedImgCh3LT)
%     title([name ' Channel 1 lifetime'],'Interpreter','none')
    colormap(jet);
    caxis(gca,scale3);
    h0 = colorbar;
    % %     ylabel(h0, ['Lifetime CH', int2str(dest_channel),' (ns)'])
    h0.Label.String = 'Lifetime (ns)';
    set(gca,'FontSize',12)
    set(gca,'LooseInset',get(gca,'TightInset'))
    frame = getframe(gcf);
    writeVideo(w,frame);
    
end
overlayTemp = overlayCh3All{end};
BW = sum(overlayTemp,3);
mask= cat(3,BW,BW,BW);
augmentedFMImgCh3LT = FM;
augmentedFMImgCh3LT(~(mask == 0)) = alpha*overlayTemp(~(mask == 0) ) + (1-alpha)*augmentedFMImgCh3LT( ~(mask == 0) );
augmentedFMImgCh3LT = augmentedFMImgCh3LT(ROI_x,ROI_y,:);
% % show image
imshow(augmentedFMImgCh3LT)
%     title([name ' Channel 1 lifetime'],'Interpreter','none')
colormap(jet);
caxis(gca,scale3);
h0 = colorbar;
% %     ylabel(h0, ['Lifetime CH', int2str(dest_channel),' (ns)'])
h0.Label.String = 'Lifetime (ns)';
set(gca,'FontSize',12)
set(gca,'LooseInset',get(gca,'TightInset'))
frame = getframe(gcf);
writeVideo(w,frame);
writeVideo(w,frame);
close(w)
end
%%
cd(root)
