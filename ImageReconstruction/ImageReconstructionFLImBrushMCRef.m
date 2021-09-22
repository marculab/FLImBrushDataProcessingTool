% code to interplocation of V6 data
%%
close all
clear
clc

% addpath(genpath(pwd))
addpath(genpath('..\'))
%% load in data
% root = 'D:\LoaclData\5ALAFLImTest\Subject031_20210804';
root = 'D:\LoaclData\5ALAFLImTest\Subject031_20210804';
%% set save path
%savePath = 'D:\LoaclData\5ALAFLImTest\Subject031_20210804\_OverlayWhiteLightL';
savePath = 'D:\LoaclData\5ALAFLImTest\Subject031_20210804\ImgMotionCorrection';
[DeConfile,DeConpath] = uigetfile([root '\*.mat'],'Please select DeCon file','MultiSelect','off');


    fileTemp = DeConfile;
    [~,runName,~] = fileparts(fileTemp);
    runNum = strsplit(runName,'_');
    runNum = str2double(runNum{5});
    %% get image file
    [imgFile,imgPath] = uigetfile([root '\*.png'],'Please select image file','MultiSelect','off');
       
    %% get txt file
    [txtFile,txtPath] = uigetfile([root '\*.txt'],'Please select txt file','MultiSelect','off');
    
    %% get location file
    [locFile,locPath] = uigetfile([root '\*.mat'],'Please select location file','MultiSelect','off');
    
    %% load in data
    load(fullfile(DeConpath,fileTemp))
    load(fullfile(locPath,locFile))
    VideoData = importVideoTextFile(fullfile(txtPath,txtFile));
    RepRate = Ch1DataObj.laserRepRate;
    
    %% open video and get the last image for augmentation
    im = imread(fullfile(imgPath, imgFile));
    figure
    image(im)
    output.img = im;
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
    
    Ch1INTCorr(G1>G1Max) = 0;
    Ch1LT(G1>G1Max) = 0;
    Ch1SNR(G1>G1Max) = 0;
    
    Ch2INTCorr(G2>G2Max) = 0;
    Ch2LT(G2>G2Max) = 0;
    Ch2SNR(G2>G2Max) = 0;
    
    Ch3INTCorr(G3>G3Max) = 0;
    Ch3LT(G3>G3Max) = 0;
    Ch3SNR(G3>G3Max) = 0;
    
    %% repopulate output data mat with interplation
    % reset start time to 0
    VideoData(:,9) = VideoData(:,9)-VideoData(1,9);
    frameTime = VideoData(:,9);
    % use CNN segmentation location data
    try
        VideoData(:,6:7) = posMC{end};
    catch
        VideoData(:,6:7) = double(flimDataset(:,7:8));
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
    [filepath,name,ext] = fileparts(fileTemp);
    save([name '_ImgReconMCRef.mat'],'output')
    disp('Reconstructed image .mat file saved successfully!')
    close all
    %% set position daya
    posData.px = xx;
    posData.py = yy;
    posData.frames = frameIdx;
    %%
    cd(savePath)
    % replotVideo(['videos\' videoName '.avi'], [name '_interp.mat'])
    radius = 7.5;
    alpha = 0.5;
    ROI_x = 1:size(im,1);
    ROI_y = 1:size(im,2);
    
    %-------------------------------------------Channel 1 lifetime----------------------------------------------------------
    figure('units','normalized','outerposition',[0.1 0.1 0.8 0.8])
    tiledlayout(2,3)
    nexttile
    scale = [floor(quantile(Ch1LT,0.1)) ceil(quantile(Ch1LT,0.9))];
    if scale(1) == scale(2)
        scale(2) = scale(1)+1;
    end
    % scale = [1 7];
    [augmentedImgCh1LT,~] = AugmentImg(im, posData, Ch1LT, scale, radius, alpha);
    augmentedImgCh1LT = augmentedImgCh1LT(ROI_x,ROI_y,:);
    % % show image
    imshow(augmentedImgCh1LT)
    title([name ' Channel 1 lifetime'],'Interpreter','none')
    colormap(jet);
    caxis(gca,scale);
    h0 = colorbar;
    % %     ylabel(h0, ['Lifetime CH', int2str(dest_channel),' (ns)'])
    h0.Label.String = 'Lifetime (ns)';
    set(gca,'FontSize',12)
    set(gca,'LooseInset',get(gca,'TightInset'))
    % exportgraphics(gca, [name '_ch1 lifetime','.jpg'],'Resolution',600);
    %----------------------------------------Channel 2 lifetime----------------------------------------------
    % figure('units','normalized','outerposition',[0 0 1 1])
    % tiledlayout(1,2)
    nexttile
    scale = [floor(quantile(Ch2LT,0.1)) ceil(quantile(Ch2LT,0.9))];
    % scale = [2 9];
    [augmentedImgCh2LT,~] = AugmentImg(im, posData, Ch2LT, scale, radius, alpha);
    augmentedImgCh2LT = augmentedImgCh2LT(ROI_x,ROI_y,:);
    % show image
    imshow(augmentedImgCh2LT)
    title([name ' Channel 2 lifetime'],'Interpreter','none')
    colormap(jet);
    caxis(gca,scale);
    h0 = colorbar;
    %     ylabel(h0, ['Lifetime CH', int2str(dest_channel),' (ns)'])
    h0.Label.String = 'Lifetime (ns)';
    set(gca,'FontSize',12)
    set(gca,'LooseInset',get(gca,'TightInset'))
    % exportgraphics(gca, [name '_ch2 lifetime','.jpg'],'Resolution',600);
    
    %------------------------------------Channel 3 lifetime------------------------------------------------------
    scale = [floor(quantile(Ch3LT,0.1)) ceil(quantile(Ch3LT,0.9))];
    % scale = [2 16];
    [augmentedImgCh3LT,~] = AugmentImg(im, posData, Ch3LT, scale, radius, alpha);
    augmentedImgCh3LT = augmentedImgCh3LT(ROI_x,ROI_y,:);
    % show image
    % figure('units','normalized','outerposition',[0 0 1 1])
    % tiledlayout(1,2)
    nexttile
    imshow(augmentedImgCh3LT)
    title([name ' Channel 3 lifetime'],'Interpreter','none')
    colormap(jet);
    caxis(gca,scale);
    h0 = colorbar;
    %     ylabel(h0, ['Lifetime CH', int2str(dest_channel),' (ns)'])
    h0.Label.String = 'Lifetime (ns)';
    set(gca,'FontSize',12)
    set(gca,'LooseInset',get(gca,'TightInset'))
    % exportgraphics(gca, [name '_ch3 lifetime','.jpg'],'Resolution',600);
    
    
    % %-------------------------------------------Channel 1 Gain----------------------------------------------------------
    scale  = [floor(quantile(G1,0.10)) ceil(quantile(G1,0.90))];
    [augmentedImgCh1G,~] = AugmentImg(im, posData, G1, scale, radius, alpha);
    augmentedImgCh1G = augmentedImgCh1G(ROI_x,ROI_y,:);
    % show image
    nexttile
    imshow(augmentedImgCh1G)
    title([name ' Channel 1 Gain'],'Interpreter','none')
    colormap(jet);
    caxis(gca,scale);
    h0 = colorbar;
    %     ylabel(h0, ['Lifetime CH', int2str(dest_channel),' (ns)'])
    h0.Label.String = 'Gain (a.u.)';
    set(gca,'FontSize',12)
    set(gca,'LooseInset',get(gca,'TightInset'))
    % exportgraphics(gca, [name '_ch1 gain','.jpg'],'Resolution',600);
    % exportgraphics(gcf, [name '_ch1','.jpg'],'Resolution',600);
    %-------------------------------------------Channel 2 Gain----------------------------------------------------------
    % scale = [mean(G2)-1*std(G2) mean(G2)+1*std(G2)];
    scale  = [floor(quantile(G2,0.10)) ceil(quantile(G2,0.90))];
    % scale(scale<0) = 0;
    [augmentedImgCh2G,~] = AugmentImg(im, posData, G2, scale, radius, alpha);
    augmentedImgCh2G = augmentedImgCh2G(ROI_x,ROI_y,:);
    % show image
    nexttile
    imshow(augmentedImgCh2G)
    title([name ' Channel 2 Gain'],'Interpreter','none')
    colormap(jet);
    caxis(gca,scale);
    h0 = colorbar;
    %     ylabel(h0, ['Lifetime CH', int2str(dest_channel),' (ns)'])
    h0.Label.String = 'Gain (a.u.)';
    set(gca,'FontSize',12)
    set(gca,'LooseInset',get(gca,'TightInset'))
    % exportgraphics(gcf, [name '_ch2','.jpg'],'Resolution',600);
    
    
    %-------------------------------------------Channel 3 Gain----------------------------------------------------------
    scale  = [floor(quantile(G3,0.10)) ceil( quantile(G3,0.90))];
    [augmentedImgCh3G,~] = AugmentImg(im, posData, G3, scale, radius, alpha);
    augmentedImgCh3G = augmentedImgCh3G(ROI_x,ROI_y,:);
    % show image
    nexttile
    imshow(augmentedImgCh3G)
    title([name ' Channel 3 Gain'],'Interpreter','none')
    colormap(jet);
    caxis(gca,scale);
    h0 = colorbar;
    %     ylabel(h0, ['Lifetime CH', int2str(dest_channel),' (ns)'])
    h0.Label.String = 'Gain (a.u.)';
    set(gca,'FontSize',12)
    set(gca,'LooseInset',get(gca,'TightInset'))
    exportgraphics(gcf, [name,'MCRef.tif'],'Resolution',600);
    
    % close all
    
    cd(root)

disp('All DONE')