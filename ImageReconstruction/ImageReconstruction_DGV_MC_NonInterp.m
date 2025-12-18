% code to interplocation of V6 data
%%
close all
clear all
clc

addpath(genpath(pwd))

%% load in data
root = 'C:\Users\Xiangnan\Box\Prostate Pilot\RARPDatabase\RARP_P3';
[DeConfile,DeConpath] = uigetfile([root '\*.mat'],'Please select DeCon file');
[Txtfile,Txtpath] = uigetfile([root '\*.txt'],'Please select text file');
[Posfile,Pospath] = uigetfile([root '\*.xlsx'],'Please select MC Pos file');
[MCfile,MCpath] = uigetfile([root '\*.csv'],'Please select MC Pos file');
[videoName,videoPath] = uigetfile([root '\*.avi'],'Please select video file');
savePath = 'C:\Users\Xiangnan\Box\Prostate Pilot\RARPDatabase\RARP_P3\AugmentedImages';
%% load in data
load(fullfile(DeConpath,DeConfile))
[ABFrame, ~, ~] = importPos(fullfile(Pospath,Posfile));
[Frame, PointIndex, X, Y] = importAllFrame(fullfile(MCpath,MCfile));
VideoData = importVideoTextFile(fullfile(Txtpath,Txtfile));
RepRate = Ch1DataObj.laserRepRate;
ABFrame = ABFrame+1;
Frame = Frame+1;
%% set up reference frame and extract location data
% refFrame = 1384;
refFrame = 1332;
PointIndex = PointIndex(Frame==refFrame);
ABFrameIdx = ABFrame(PointIndex+1);
X = X(Frame==refFrame);
Y = Y(Frame==refFrame);


%% open video and get the last image for augmentation
v = VideoReader(fullfile(videoPath, videoName));
im = read(v,refFrame);
imx = size(im,2);
imy = size(im,1);
im = rgb2gray(im);
im = repmat(im,[1 1 3]);
% im = imresize(im, [720 1280]);
% im = read(v,1);
figure
imshow(im)
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

Ch4INTCorr = circshift(Ch4DataObj.Lg_INTsGainCorrected,shift);
Ch4LT = circshift(Ch4DataObj.Lg_LTs,shift);
Ch4SNR = circshift(Ch4DataObj.SNR,shift)';
G4 = circshift(Ch4DataObj.gain,shift);

%% get FLIm data
FLImIdx = VideoData(ABFrameIdx,1);

%% save data
% cd(DeConpath)
[filepath,name,ext] = fileparts(videoName);
% save([name '_ImgRecon.mat'],'output')
% disp('Reconstructed image .mat file saved successfully!')
% close all
%% set position daya

%% filter data
position = [X Y];
Ch1LT = Ch1LT(FLImIdx);
Ch2LT = Ch2LT(FLImIdx);
Ch3LT = Ch3LT(FLImIdx);
Ch4LT = Ch4LT(FLImIdx);
G1 = G1(FLImIdx);
G2 = G2(FLImIdx);
G3 = G3(FLImIdx);
G4 = G4(FLImIdx);

gainMask = (G1<200)&(G2<1000)&(G3<1000)&(position(:,1)>386);
% gainMask(1950:end) = [];
position = position(gainMask,:);
Ch1LT = Ch1LT(gainMask);
Ch2LT = Ch2LT(gainMask);
Ch3LT = Ch3LT(gainMask);
Ch4LT = Ch4LT(gainMask);
G1 = G1(gainMask);
G2 = G2(gainMask);
G3 = G3(gainMask);
G4 = G4(gainMask);


%%
cd(savePath)
% replotVideo(['videos\' videoName '.avi'], [name '_interp.mat'])
radius = 25;
alpha = 0.5;
cmap = jet(256);
sizeArr = ones(size(Ch1LT))*radius;

%-------------------------------------------Channel 1 lifetime----------------------------------------------------------
figure('units','normalized','outerposition',[0 0 1 1])
tiledlayout(1,2)
nexttile
scale = [floor(quantile(Ch1LT,0.1)) ceil(quantile(Ch1LT,0.9))];
% scale = [2 5];
% [augmentedImg,~] = AugmentImg(im, posData, Ch1LT, scale, radius, alpha, cmap);
augmentedImg = Overlay_DGV(size(im,2), size(im,1), sizeArr, position, Ch1LT, im, Ch1SNR, 0, scale);

% show image
imshow(augmentedImg); 
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
augmentedImg = Overlay_DGV(size(im,2), size(im,1), sizeArr, position, G1, im, Ch1SNR, 0,scale);
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
exportgraphics(gcf, [name '_ch1_DGV','.jpg'],'Resolution',600);
%----------------------------------------Channel 2 lifetime----------------------------------------------
figure('units','normalized','outerposition',[0 0 1 1])
tiledlayout(1,2)
nexttile

scale = [floor(quantile(Ch2LT,0.1)) ceil(quantile(Ch2LT,0.9))];
% scale = [2 5];
augmentedImg = Overlay_DGV(size(im,2), size(im,1), sizeArr, position, Ch2LT, im, Ch1SNR, 0, scale);
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
augmentedImg = Overlay_DGV(size(im,2), size(im,1), sizeArr, position, G2, im, Ch1SNR, 0, scale);
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
exportgraphics(gcf, [name '_ch2_DGV','.jpg'],'Resolution',600);

%------------------------------------Channel 3 lifetime------------------------------------------------------
scale = [floor(quantile(Ch3LT,0.1)) ceil(quantile(Ch3LT,0.9))];
% scale = [2 5];
augmentedImg = Overlay_DGV(size(im,2), size(im,1), sizeArr, position, Ch3LT, im, Ch1SNR, 0, scale);
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
augmentedImg = Overlay_DGV(size(im,2), size(im,1), sizeArr, position, G3, im, Ch1SNR, 0, scale);
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
exportgraphics(gcf, [name '_ch3_DGV','.jpg'],'Resolution',600);

%------------------------------------Channel 4 lifetime------------------------------------------------------
scale = [floor(quantile(Ch4LT,0.1)) ceil(quantile(Ch4LT,0.9))];
% scale = [2 5];
augmentedImg = Overlay_DGV(size(im,2), size(im,1), sizeArr, position, Ch4LT, im, Ch1SNR, 0, scale);
% show image
figure('units','normalized','outerposition',[0 0 1 1])
tiledlayout(1,2)
nexttile
imshow(augmentedImg)
title([name ' Channel 4 lifetime'],'Interpreter','none')
colormap(jet);
caxis(gca,scale);
h0 = colorbar;
%     ylabel(h0, ['Lifetime CH', int2str(dest_channel),' (ns)'])
h0.Label.String = 'Lifetime (ns)';
set(gca,'FontSize',15)
set(gca,'LooseInset',get(gca,'TightInset'))
% exportgraphics(gca, [name '_ch3 lifetime','.jpg'],'Resolution',600);

%-------------------------------------------Channel 4 Gain----------------------------------------------------------
scale  = [floor(quantile(G4,0.10)) ceil( quantile(G4,0.90))];
augmentedImg = Overlay_DGV(size(im,2), size(im,1), sizeArr, position, G4, im, Ch1SNR, 0, scale);
% show image
nexttile
imshow(augmentedImg)
title([name ' Channel 4 Gain'],'Interpreter','none')
colormap(jet);
caxis(gca,scale);
h0 = colorbar;
%     ylabel(h0, ['Lifetime CH', int2str(dest_channel),' (ns)'])
h0.Label.String = 'Gain (a.u.)';
set(gca,'FontSize',15)
set(gca,'LooseInset',get(gca,'TightInset'))
exportgraphics(gcf, [name '_ch4_DGV','.jpg'],'Resolution',600);

% close all