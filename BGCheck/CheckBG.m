% code to check FLImBRUSH background file, Xiangnan, 05-12-2022
%%
clear 
close all
clc

%% add BG class to path
addpath(genpath('..\Algorithms'))
%% load in APD detector file
% channel 1
[APD1file,APD1path] = uigetfile('.mat','Please select APD detector file','MultiSelect',"off");
APD1 = load(fullfile(APD1path,APD1file));
% channel 2
[APD2file,APD2path] = uigetfile('.mat','Please select APD detector file','MultiSelect',"off");
APD2 = load(fullfile(APD2path,APD2file));
% channel 3
[APD3file,APD3path] = uigetfile('.mat','Please select APD detector file','MultiSelect',"off");
APD3 = load(fullfile(APD3path,APD3file));

%% select background file
[BGfile,BGpath] = uigetfile('.tdms','Please select BG tdms file','MultiSelect',"off");
BGObj = backGround(fullfile(BGpath,BGfile));
BGObj.loadBG;
BG1_Gain = interp1(APD1.gainV,APD1.gain,BGObj.CtrlV1);
BG2_Gain = interp1(APD2.gainV,APD2.gain,BGObj.CtrlV2);
BG3_Gain = interp1(APD3.gainV,APD3.gain,BGObj.CtrlV3);

figure
plot(BGObj.bgCh1)
title(['Channel 1 Background, Gain = ' num2str(BG1_Gain,'%.2f')]);

figure
plot(BGObj.bgCh2)
title(['Channel 2 Background, Gain = ' num2str(BG2_Gain,'%.2f')]);

figure
plot(BGObj.bgCh3)
title(['Channel 3 Background, Gain = ' num2str(BG3_Gain,'%.2f')]);
