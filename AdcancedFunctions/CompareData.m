clear
close all
clc

addpath(genpath('E:\MyGitRepo\FLImBrushDataProcessingTool\Algorithms'))
%%
% path1 = ['Z:\Da Vinci Robot Study (100 patients) no iRF cap no CFDBug' ...
%     '\Data_100_Patient_Study\Subject_146\Triplex\ALL_DECONVOLVED_FILES' ...
%     '\P800853_07_15_22_02_12.5GS_IRF_lite.mat'];
% 
% path2 = ['Z:\Da Vinci Robot Study (100 patients) no iRFcap wCFDBug' ...
%     '\Data_100_Patient_Study\Subject_146\Triplex\ALL_DECONVOLVED_FILES' ...
%     '\P800853_07_15_22_02_12.5GS_wCFDBug_lite.mat'];

path1 = 'E:\Patient175\CorrectCFD\P800853_03_29_23_02_12.5GS_IRF';
path2 = 'E:\Patient175\WrongCFD\P800853_03_29_23_02_12.5GS_wCFDBug';
% path1 = 'E:\Patient146\CorrectCFD\P800853_07_15_22_02_12.5GS_IRF';
% path2 = 'E:\Patient146\WrongCFD\P800853_07_15_22_02_12.5GS_wCFDBug';

% path1 = 'E:\DataProcessingRelated\20221108FB_vs_V4\20221108FBvsV4\Decon 2024-7-17 13-11\20221108FBvsV4_01.mat';
% path2 = 'E:\DataProcessingRelated\20221108FB_vs_V4\20221108FBvsV4\Decon 2024-7-17 13-13_SingleShift\20221108FBvsV4_01.mat';

% path1 = 'E:\C440\CorrectCFD\20221108FBvsV4_01';
% path2 = 'E:\C440\WrongCFD\20221108FBvsV4_01';
%% load in data
data1 = load(path1);
data2 = load(path1);

%% check truncated data
dataT1 = data1.Ch2DataObj.dataT;
dataT2 = data2.Ch2DataObj.dataT;
idx = 450:455;
figure
tiledlayout(1,2)
nexttile
plot(dataT1(:,idx))
nexttile
plot(dataT2(:,idx))
%% redo CFD align dataT
dataT1 = alignWaveform_CFDNew(dataT1, 2.4, 0.08,0.5);
dataT2 = alignWaveform_CFDNew(dataT2, 2.4, 0.08,0.5);
figure
tiledlayout(1,2)
nexttile
plot(dataT1(:,idx))
nexttile
plot(dataT2(:,idx))
%% rerun deonvolution
% data1 = load('E:\DataProcessingRelated\20221108FB_vs_V4\20221108FBvsV4\Decon 2024-3-14 10-36 fixed CFD\20221108FBvsV4_01.mat')

runDeconLG(data1.Ch1DataObj,[450:550],12,0.916); 
runDeconLG(data2.Ch1DataObj,[450:550],12,0.916); 
runDeconLG(data1.Ch2DataObj,[450:550],12,0.916); 
runDeconLG(data2.Ch2DataObj,[450:550],12,0.916); 
runDeconLG(data1.Ch3DataObj,[450:550],12,0.916); 
runDeconLG(data2.Ch3DataObj,[450:550],12,0.916); 
%%
figure
tiledlayout(1,3)
nexttile
scatter(data1.Ch1DataObj.Lg_LTs,data2.Ch1DataObj.Lg_LTs,'b.')
hold on
plot([0 10],[0 10],'r--')
xlim([3 6])
ylim([3 6])
title('Channel 1 LT')

nexttile
scatter(data1.Ch2DataObj.Lg_LTs,data2.Ch2DataObj.Lg_LTs,'b.')
hold on
plot([0 10],[0 10],'r--')
% xlim([0 8])
% ylim([0 8])
% axis equal
xlim([3 6])
ylim([3 6])
xlabel('Correct CFD')
ylabel('Wrong CFD')
title('Channel 2 LT')

nexttile
scatter(data1.Ch3DataObj.Lg_LTs,data2.Ch3DataObj.Lg_LTs,'b.')
hold on
plot([0 10],[0 10],'r--')
xlim([3 6])
ylim([3 6])
title('Channel 3 LT')

%% plot histogram from both run
figure
tiledlayout(1,3)
bin = 2:0.1:6;
nexttile
histogram(data1.Ch1DataObj.Lg_LTs,bin,'FaceColor','r')
hold on
histogram(data2.Ch1DataObj.Lg_LTs,bin,'FaceColor','g')
hold off
% xlim([0 8])
% ylim([0 8])
% title('Channel 1 LT')

nexttile
histogram(data1.Ch2DataObj.Lg_LTs,bin,'FaceColor','r')
hold on
histogram(data2.Ch2DataObj.Lg_LTs,bin,'FaceColor','g')
hold off

nexttile
histogram(data1.Ch3DataObj.Lg_LTs,bin,'FaceColor','r')
hold on
histogram(data2.Ch3DataObj.Lg_LTs,bin,'FaceColor','g')
hold off
% xlim([0 8])
% ylim([0 8])
% title('Channel 3 LT')
%% plot histogram difference
figure
tiledlayout(1,3)
nexttile
histogram(data1.Ch1DataObj.Lg_LTs-data2.Ch1DataObj.Lg_LTs,[-0.2:0.005:0.2])
% xlim([0 8])
% ylim([0 8])
% title('Channel 1 LT')

nexttile
histogram(data1.Ch2DataObj.Lg_LTs-data2.Ch2DataObj.Lg_LTs,[-0.2:0.005:0.2])
% xlim([0 8])
% ylim([0 8])
% title('Channel 2 LT')

nexttile
histogram(data1.Ch3DataObj.Lg_LTs-data2.Ch3DataObj.Lg_LTs,[-0.2:0.005:0.2])
% xlim([0 8])
% ylim([0 8])
% title('Channel 3 LT')

%% plot histogram
figure
tiledlayout(1,3)
nexttile
scatter(data1.Ch1DataObj.Lg_LTs,data1.Ch1DataObj.Lg_LTs-data2.Ch1DataObj.Lg_LTs,'b.')
xlim([0 8])
% ylim([0 8])
% title('Channel 1 LT')

nexttile
scatter(data1.Ch2DataObj.Lg_LTs,data1.Ch2DataObj.Lg_LTs-data2.Ch2DataObj.Lg_LTs,'b.')
xlim([0 8])
% ylim([0 8])
% title('Channel 2 LT')

nexttile
scatter(data1.Ch3DataObj.Lg_LTs,data1.Ch3DataObj.Lg_LTs-data2.Ch3DataObj.Lg_LTs,'b.')
xlim([0 8])
% ylim([0 8])
% title('Channel 3 LT')
%% plot lifetime difference trace

Lg_LT_diff = data1.Ch2DataObj.Lg_LTs-data2.Ch2DataObj.Lg_LTs;
figure
plot(Lg_LT_diff)

figure
scatter(Lg_LT_diff,data1.Ch2DataObj.irfIdx)


%% get aligned waveform and fit
data1Aligned = get(data1.Ch2DataObj,'wf_aligned');
data2Aligned = get(data2.Ch2DataObj,'wf_aligned');
data1Fit = get(data1.Ch2DataObj,'fit');
data2Fit = get(data2.Ch2DataObj,'fit');
decay1 = get(data1.Ch2DataObj,'decay');
decay2 = get(data2.Ch2DataObj,'decay');
%% correlation between waveform max and lifetime difference
[~,maxIdx1] = max(data1Aligned);
[~,maxIdx2] = max(data2Aligned);
maxDiff = maxIdx1-maxIdx2;
figure
scatter(maxDiff,Lg_LT_diff)
xlim([-1 2])
intersect(find(maxDiff==0),find(Lg_LT_diff<-0.2))
%% get iRF
irf1 = data1.Ch2DataObj.APDObj.irfTNorm(:,data1.Ch2DataObj.irfIdx);
irf2 = data2.Ch2DataObj.APDObj.irfTNorm(:,data2.Ch2DataObj.irfIdx);
%---------------compute AMD lifetime---------------------------------------
dt = 0.08;
t = 1:1:size(data1.Ch2DataObj.dataT,1);
t = t*dt;
AMDIRF1 = irf1(:,:)'*t'./sum(irf1(:,:))';
% AMDIRF2 = irf2(:,idx)'*t'/sum(irf2(:,idx))
% AMDLT1 = data1.Ch2DataObj.dataT(:,:)'*t'./sum(data1.Ch2DataObj.dataT(:,:))'-AMDIRF1
% AMDLT2 = data2.Ch2DataObj.dataT(:,:)'*t'./sum(data2.Ch2DataObj.dataT(:,:))'-AMDIRF1
% figure
% histogram(AMDLT1,[2:0.05:8])
% 
% figure
% histogram(data1.Ch2DataObj.Lg_LTs,[2:0.05:8])

AMDLT1 = data1Aligned(:,:)'*t'./sum(data1Aligned(:,:))'-AMDIRF1;
AMDLT2 = data2Aligned(:,:)'*t'./sum(data2Aligned(:,:))'-AMDIRF1;
%% look for specific point
% idx = 7009;
idx = 5008;
%-----------------------plot truncated waveform----------------------------
figure
plot(data1.Ch2DataObj.dataT(:,idx))
hold on
plot(data2.Ch2DataObj.dataT(:,idx))
hold off
%----------------------------plot fitting ----------------------------
figure
plot(data1Aligned(:,idx))
hold on
plot(data1Fit(:,idx))
hold off

figure
plot(data2Aligned(:,idx))
hold on
plot(data2Fit(:,idx))
hold off
%-------------------------plot irf aligned waveform------------------------
figure
tiledlayout(1,3)
nexttile
plot(data1Aligned(:,idx))
hold on
plot(data2Aligned(:,idx))
plot(irf1(:,idx)/max(irf1(:,idx))*0.5)
plot(irf2(:,idx)/max(irf2(:,idx))*0.5)
hold off
grid on
nexttile
plot(decay1(:,idx))
hold on
plot(decay2(:,idx))
hold off
grid on
title('Decay')
nexttile
plot(decay1(:,idx)/max(decay1(:,idx)))
hold on
plot(decay2(:,idx)/max(decay2(:,idx)))
hold off
grid on
title('Normalized decay')
data1.Ch2DataObj.shift(idx)
data2.Ch2DataObj.shift(idx)
data1.Ch2DataObj.Lg_LTs(idx)
data2.Ch2DataObj.Lg_LTs(idx)
AMDLT1(idx)
AMDLT2(idx)
%%

histogram(AMDLT1,[2:0.05:8])
figure
scatter(AMDLT1,AMDLT2,'b.')
xlim([2 5])
ylim([2 5])
title('Channel 2 LT')
hold on
plot([0 10],[0 10],'r-.')
hold off

figure
histogram(AMDLT1-AMDLT2,[-0.5:0.05:0.5])


data1.Ch2DataObj.Lg_LTs(idx)
data2.Ch2DataObj.Lg_LTs(idx)
AMDLT1(idx)
AMDLT2(idx)
%--------------------plot Laguerre vs AMD lifetime------------------------
figure
tiledlayout(1,2)
nexttile
scatter(data1.Ch2DataObj.Lg_LTs,AMDLT1,'b.')
hold on
plot([0 10],[0 10],'r-.')
hold off
axis equal
title('Correct CFD')
xlim([1 7])
ylim([1 7])
nexttile
scatter(data2.Ch2DataObj.Lg_LTs,AMDLT1,'b.')
hold on
plot([0 10],[0 10],'r-.')
hold off
axis equal
title('Wrong CFD')
xlim([1 7])
ylim([1 7])
% data1Aligned = get(data1.Ch2DataObj,'wf_aligned');

% figure
% plot(data1Aligned(:,idx))
% hold on
% plot(circshift(data2Aligned(:,idx),1))
% hold off
% AMDLT1 = data1Aligned(:,:)'*t'/sum(data1Aligned(:,:))'
% AMDLT2 = circshift(data2Aligned(:,idx),0)'*t'/sum(circshift(data2Aligned(:,idx),0))

% data2Aligned = get(data2.Ch2DataObj,'wf_aligned');

figure
plot(data1Fit(:,idx))
hold on
plot(data2Fit(:,idx))
hold off

figure
plot(decay1(:,idx))
hold on
plot(decay2(:,idx))
hold off
title('Decay')
sum(decay1(:,idx))
sum(decay2(:,idx))

figure
plot(decay1(:,idx)/max(decay1(:,idx)))
hold on
plot(decay2(:,idx)/max(decay2(:,idx)))
hold off
title('Normalized decay')

ROI = 1:680;
AMDLT1 = data1Fit(ROI,idx)'*t(ROI)'/sum(data1Fit(ROI,idx))
AMDLT2 = data2Fit(ROI,idx)'*t(ROI)'/sum(data2Fit(ROI,idx))


data1.Ch2DataObj.stat_test.chi2.stat(idx)
data2.Ch2DataObj.stat_test.chi2.stat(idx)

figure
histogram(data1.Ch2DataObj.stat_test.chi2.stat,[0:0.05:6],'FaceColor','r')
hold on
histogram(data2.Ch2DataObj.stat_test.chi2.stat,[0:0.05:6],'FaceColor','g')
hold off

figure
scatter(data1.Ch2DataObj.stat_test.chi2.stat,data2.Ch2DataObj.stat_test.chi2.stat,'b.')
hold on
plot([0,10],[0 10],'r--')

%% 
idx = 4009;
figure
plot([data1Aligned(:,idx) data2Aligned(:,idx)])
% idx = 4909;
dataNew = [data1.Ch2DataObj.dataT(:,idx) data2.Ch2DataObj.dataT(:,idx)];
irfNew = irf1(:,idx);
channelDataStruct = ChannelData(dataNew,irfNew,0.08,[],[],[],[]);
laguerreObj = LaguerreModel(channelDataStruct,12, 0.916);
laguerreObj.estimate_laguerre([450:550],-5:20);
laguerreObj.LTs
laguerreObj.shift
plot(laguerreObj.spec_aligned)