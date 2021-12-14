% code to run multiexpoential fit to FLImBRUSH output

%%
clear
close all
clc
% msgbox({'V4 processing code is needed!'; 'Please update the path to V4 processing code in section 2 of the code!'})
%% add V4 processing code to path, please update the path accordingly
% addpath(genpath('./MultiExp Generic')) %add path to exponential code
addpath(genpath('C:\Users\Xiangnan\Documents\MyGitRepo\TRFS_TRIPLEX_GUI\Algorithms'))
addpath(genpath('C:\Users\Xiangnan\Documents\MyGitRepo\TRFS_TRIPLEX_GUI\Preprocessing'))


%% select file
root = 'D:\LoaclData\MarginsForBinning';
[DeConfile,DeConpath] = uigetfile([root '\*.mat'],'Please select DeCon file','MultiSelect','on');
if iscell(DeConfile)
    numOfFiles = length(DeConfile);
else
    numOfFiles = 1;
end
for k = 1:numOfFiles
    % load in file
    if numOfFiles==1
        P = fullfile(DeConpath,DeConfile);
        [filepath,name,ext]=fileparts(DeConfile);
    else
        P = fullfile(DeConpath,DeConfile{k});
        [filepath,name,ext]=fileparts(DeConfile{k});
    end
    load(P)
    name
    %% SNR filter data
    dataObj = Data{2}; %% set up your channel here
    WFAligned = dataObj.channeldata.data;
    SNR = dataObj.channeldata.SNR;
    SNRIdx = SNR<30;
    WFAligned(:,SNRIdx) = [];
    irf = dataObj.channeldata.iIRF;
    numOfWF = size(WFAligned,2);
    n_max = floor(log2(numOfWF));
    WFAligned = WFAligned(:,end-2^n_max+1:end); % pick 2^n_max waveforms
    %%
    %     figure('units','pixels','outerposition',[100 100 1920 1080])
    %     tiledlayout(3,3)
    
    for c = n_max
        % get Laguerre result
        wfAvg = power(2,c);
        numOfAveWf = size(WFAligned,2)/wfAvg;
        WF = zeros(size(WFAligned,1),numOfAveWf);
        
        for i = 1:numOfAveWf
            WF(:,i) = mean(WFAligned(:,wfAvg*(i-1)+1:wfAvg*i),2);
        end
        
        SNR_new = 20*log10(max(WF)./std(WF(1:40,:)));
        %% bi-exponential fit
        % order2 = 2;
        % [A2, T2, avglife2, intensity2, fit2, raw2, decay2] = multiexp_fit(WF,0.08,irf,order2,[]);
        % res2 = WF-fit2;
        % res2SE = sum(res2.^2);
        % LT2Decay = h_lifet(decay2,0.08);
        % W2 = A2./sum(A2);
        % WTemp = WTemp';
        
        %% tri-exponential fit
        order3 = 3;
        [A3, T3, avglife3, intensity3, fit3, raw3, decay3] = multiexp_fit(WF,0.08,irf,order3,[]);
        res3 = WF-fit3;
        res3SE = sum(res3.^2);
        LT3Decay = h_lifet(decay3,0.08);
        W3 = A3./sum(A3);
        FC3 = A3.*T3./sum(A3.*T3);
        
        
    end
    dataForPhasor.rawWF = WFAligned;
    dataForPhasor.iRF = irf;
    dataForPhasor.a = A3;
    dataForPhasor.tau = T3;
    dataForPhasor.frac_Contribution = FC3;
    
    save(fullfile(DeConpath,[name,'_Ch2_Phasor.mat']),'dataForPhasor');
    
    
    %     exportgraphics(gcf,fullfile(DeConpath,[name,'Ch2.tiff']),'Resolution',600)
end

close all
