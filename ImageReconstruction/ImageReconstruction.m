% code to interplocation of V6 data
%%
close all
clear all
clc

addpath(genpath(pwd))

%% load in data
root = '';
[DeConfile,DeConpath] = uigetfile([root '\*.mat'],'Please select DeCon file');
[Txtfile,Txtpath] = uigetfile([root '\*.txt'],'Please select text file');
[Posfile,Pospath] = uigetfile([root '\*.mat'],'Please select CNN Segementation Pos file');
[videoName,videoPath] = uigetfile([root '\*.avi'],'Please select video file');
savePath = 'C:\Users\Xiangnan\Box Sync\FLImBrush vs V4\20210505\FLImBrush\20210505_neuro024\AugmentedImages';
%% load in data
load(fullfile(DeConpath,DeConfile))
load(fullfile(Pospath,Posfile))
VideoData = importVideoTextFile(fullfile(Txtpath,Txtfile));
RepRate = Ch1DataObj.laserRepRate;

%% open video and get the last image for augmentation
v = VideoReader(fullfile(videoPath, videoName));
im = read(v,inf);
% im = read(v,1);
figure
image(im)
output.img = im;
%% calculate num of data points
% shift = length(Ch1LT)-(MetaData(end,9)-MetaData(1,9))/1000*120/4;
% shift = round(shift);
shift=0; %shift need to be 0!

Ch1INTCorr = circshift(Ch1DataObj.Lg_INTsGainCorrected,shift);
Ch1LT = circshift(Ch1DataObj.Lg_LTs,shift);
Ch1SNR = circshift(Ch1DataObj.SNR,shift)';
G1 = circshift(Ch1DataObj.gain,shift);

Ch2INTCorr = circshift(Ch2DataObj.Lg_INTsGainCorrected,shift);
Ch2LT = circshift(Ch2DataObj.Lg_LTs,shift);
Ch2SNR = circshift(Ch2DataObj.SNR,shift)';
G2 = circshift(Ch2DataObj.gain,shift);

Ch3INTCorr = circshift(Ch3DataObj.Lg_INTsGainCorrected,shift);
Ch3LT = circshift(Ch3DataObj.Lg_LTs,shift);
Ch3SNR = circshift(Ch3DataObj.SNR,shift)';
G3 = circshift(Ch3DataObj.gain,shift);

%% repopulate output data mat with interplation
% reset start time to 0
VideoData(:,9) = VideoData(:,9)-VideoData(1,9);
% use CNN segmentation location data
try
    VideoData(:,6:7) = double(pos(:,7:8));
catch
    VideoData(:,6:7) = double(pos_new(:,1:2));
end
% Length = 70000;
Length = VideoData(end,9);
VideoData(VideoData(:,9)>Length,:) = [];
numOfVideoFrame = size(VideoData,1);
removeFlag = zeros(numOfVideoFrame,1);
%------------remove single 0 lines---------------------------------------%
for i = 1:numOfVideoFrame
    posCurrent = VideoData(i,6)+ VideoData(i,7);
    if posCurrent==0 % if current pos is all 0 check before and after
        if i==1
            posBefore = 1;
            posAfter = VideoData(i+1,6)+ VideoData(i+1,7);
        else if i==numOfVideoFrame
                posBefore = VideoData(i-1,6)+ VideoData(i-1,7);
                posAfter = 1;
            else
                posBefore = VideoData(i-1,6)+ VideoData(i-1,7);
                posAfter = VideoData(i+1,6)+ VideoData(i+1,7);
            end
            removeFlag(i) = posBefore&posAfter;
        end
    end
end
removeIdx = find(removeFlag);
VideoData(removeIdx,:) = [];

% remove all 0 points
temp = VideoData(:,6)+VideoData(:,7);
ZeroIdx = find(temp==0);
VideoData(ZeroIdx,:)=[];

frameNum = 1:size(VideoData,1);
frameNum = frameNum';

time = 0:1:size(Ch1LT,1)-1;
time = time'*1000/RepRate*4;
frameT = VideoData(:,9);
[frameT,ia] = unique(frameT); % find duplicate frame
VideoData = VideoData(ia,:); % remove duplicated data
frameNum = frameNum(ia); % remove duplicated frame
frameIdx = interp1(frameT,frameNum,time);
frameIdx = ceil(frameIdx);
output.frame = frameIdx;
VideoData = [frameNum VideoData];


%% interpolate locations
xx = interp1(VideoData(:,10),VideoData(:,7),time);
xx(isnan(xx))=0; % replace NaN with 0
yy = interp1(VideoData(:,10),VideoData(:,8),time);
yy(isnan(yy))=0; % replace NaN with 0
rr = interp1(VideoData(:,10),VideoData(:,9),time);
rr(isnan(rr)) = 0;
output.xx = xx;
output.yy = yy;
output.rr = rr;
output.lt1 = Ch1LT;
output.int1 = Ch1INTCorr;
output.lt2 = Ch2LT;
output.int2 = Ch2INTCorr;
output.lt3 = Ch3LT;
output.int3 = Ch3INTCorr;
output.snr1 = Ch1SNR;
output.snr2 = Ch2SNR;
output.snr3 = Ch3SNR;
output.gain1 = G1;
output.gain2 = G2;
output.gain3 = G3;

%% save data
cd(DeConpath)
[filepath,name,ext] = fileparts(videoName);
save([name '_ImgRecon.mat'],'output')
disp('Reconstructed image .mat file saved successfully!')
close all
%% set position daya
posData.px = xx;
posData.py = yy;
posData.frames = frameIdx;
%%
cd(savePath)
% replotVideo(['videos\' videoName '.avi'], [name '_interp.mat'])
radius = 10;
alpha = 0.5;
%-------------------------------------------Channel 1 lifetime----------------------------------------------------------
figure('units','normalized','outerposition',[0 0 1 1])
tiledlayout(1,2)
nexttile
% scale = [floor(quantile(Ch1LT,0.1)) ceil(quantile(Ch1LT,0.9))];
scale = [2 5];
[augmentedImg,~] = AugmentImg(im, posData, Ch1LT, scale, radius, alpha);
% show image
imshow(augmentedImg)
title([name ' Channel 1 lifetime'],'Interpreter','none')
colormap(jet);
caxis(gca,scale);
h0 = colorbar;
%     ylabel(h0, ['Lifetime CH', int2str(dest_channel),' (ns)'])
h0.Label.String = 'Lifetime (ns)';
set(gca,'FontSize',15)
set(gca,'LooseInset',get(gca,'TightInset'))
% exportgraphics(gca, [name '_ch1 lifetime','.jpg'],'Resolution',600);
%-------------------------------------------Channel 1 Gain----------------------------------------------------------
scale  = [floor(quantile(G1,0.10)) ceil(quantile(G1,0.90))];
[augmentedImg,~] = AugmentImg(im, posData, G1, scale, radius, alpha);
% show image
nexttile
imshow(augmentedImg)
title([name ' Channel 1 Gain'],'Interpreter','none')
colormap(jet);
caxis(gca,scale);
h0 = colorbar;
%     ylabel(h0, ['Lifetime CH', int2str(dest_channel),' (ns)'])
h0.Label.String = 'Gain (a.u.)';
set(gca,'FontSize',15)
set(gca,'LooseInset',get(gca,'TightInset'))
% exportgraphics(gca, [name '_ch1 gain','.jpg'],'Resolution',600);
exportgraphics(gcf, [name '_ch1','.jpg'],'Resolution',600);
%----------------------------------------Channel 2 lifetime----------------------------------------------
figure('units','normalized','outerposition',[0 0 1 1])
tiledlayout(1,2)
nexttile

% scale = [floor(quantile(Ch2LT,0.1)) ceil(quantile(Ch2LT,0.9))];
scale = [2 5];
[augmentedImg,~] = AugmentImg(im, posData, Ch2LT, scale, radius, alpha);
% show image
imshow(augmentedImg)
title([name ' Channel 2 lifetime'],'Interpreter','none')
colormap(jet);
caxis(gca,scale);
h0 = colorbar;
%     ylabel(h0, ['Lifetime CH', int2str(dest_channel),' (ns)'])
h0.Label.String = 'Lifetime (ns)';
set(gca,'FontSize',15)
set(gca,'LooseInset',get(gca,'TightInset'))
% exportgraphics(gca, [name '_ch2 lifetime','.jpg'],'Resolution',600);
%-------------------------------------------Channel 2 Gain----------------------------------------------------------
% scale = [mean(G2)-1*std(G2) mean(G2)+1*std(G2)];
scale  = [floor(quantile(G2,0.10)) ceil(quantile(G2,0.90))];
% scale(scale<0) = 0;
[augmentedImg,~] = AugmentImg(im, posData, G2, scale, radius, alpha);
% show image
nexttile
imshow(augmentedImg)
title([name ' Channel 2 Gain'],'Interpreter','none')
colormap(jet);
caxis(gca,scale);
h0 = colorbar;
%     ylabel(h0, ['Lifetime CH', int2str(dest_channel),' (ns)'])
h0.Label.String = 'Gain (a.u.)';
set(gca,'FontSize',15)
set(gca,'LooseInset',get(gca,'TightInset'))
exportgraphics(gcf, [name '_ch2','.jpg'],'Resolution',600);

%------------------------------------Channel 3 lifetime------------------------------------------------------
% scale = [floor(quantile(Ch3LT,0.1)) ceil(quantile(Ch3LT,0.9))];
scale = [2 5];
[augmentedImg,~] = AugmentImg(im, posData, Ch3LT, scale, radius, alpha);
% show image
figure('units','normalized','outerposition',[0 0 1 1])
tiledlayout(1,2)
nexttile
imshow(augmentedImg)
title([name ' Channel 3 lifetime'],'Interpreter','none')
colormap(jet);
caxis(gca,scale);
h0 = colorbar;
%     ylabel(h0, ['Lifetime CH', int2str(dest_channel),' (ns)'])
h0.Label.String = 'Lifetime (ns)';
set(gca,'FontSize',15)
set(gca,'LooseInset',get(gca,'TightInset'))
% exportgraphics(gca, [name '_ch3 lifetime','.jpg'],'Resolution',600);

%-------------------------------------------Channel 3 Gain----------------------------------------------------------
scale  = [floor(quantile(G3,0.10)) ceil( quantile(G3,0.90))];
[augmentedImg,~] = AugmentImg(im, posData, G3, scale, radius, alpha);
% show image
nexttile
imshow(augmentedImg)
title([name ' Channel 3 Gain'],'Interpreter','none')
colormap(jet);
caxis(gca,scale);
h0 = colorbar;
%     ylabel(h0, ['Lifetime CH', int2str(dest_channel),' (ns)'])
h0.Label.String = 'Gain (a.u.)';
set(gca,'FontSize',15)
set(gca,'LooseInset',get(gca,'TightInset'))
exportgraphics(gcf, [name '_ch3','.jpg'],'Resolution',600);

close all