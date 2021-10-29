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
    if k==1
        P = fullfile(DeConpath,DeConfile);
        [filepath,name,ext]=fileparts(DeConfile);
    else
        P = fullfile(DeConpath,DeConfile{k});
        [filepath,name,ext]=fileparts(DeConfile{k});
    end
    load(P)
    
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
    figure('units','pixels','outerposition',[100 100 1920 1080])
    tiledlayout(3,3)
    
    for c = 0%:n_max
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
        %% tri-exponential fit
        % order4 = 4;
        % [A4, T4, avglife4, intensity4, fit4, raw4, decay4] = multiexp_fit(WF,0.08,irf,order4,[]);
        % res4 = WF-fit4;
        % res4SE = sum(res4.^2);
        % LT4Decay = h_lifet(decay4,0.08);
        % W4 = A4./sum(A4);
        
        %% plot result
        nexttile
        scatter(FC3(1,:),T3(1,:),'bo','LineWidth',2)
        hold on
        scatter(FC3(2,:),T3(2,:),'g+','LineWidth',2)
        scatter(FC3(3,:),T3(3,:),SNR_new,'rs','LineWidth',2)
        hold off
        legend('exp1','exp2','exp3')
        xlabel('Fractional Contribution')
        ylabel('Lifetime (ns)')
        xlim([0 1])
        ylim([0 25])
        grid on
        grid minor
        title(sprintf('Averaging = %d',2^c));
        
        
        % idx = 1;
        %
        % % plot bi-exponential
        % nexttile
        % plot(WF(:,idx),'b.','MarkerSize',10)
        % hold on
        % plot(fit2(:,idx),'r-','LineWidth',1)
        % plot(res2(:,idx)-0.1,'m.')
        % yline(-0.1+0.025,'r--','res=0.025')
        % yline(-0.1-0.025,'r--')
        % text(150, -0.15, ['\Sigma res = ', num2str(res2SE)])
        % yticks((-0.2:0.1:1))
        % yticklabels({'','res=0','0','0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','1.0',})
        % grid on
        % grid minor
        % axis tight
        % ylim([-0.2 1])
        % legend('raw data', 'fitting', 'residue')
        % title(sprintf('Channel %d, Exponential fitting, order = %d, WF index %d,\n lifetime = %.3f(formula), lifetime = %.3f(decay),\n a1=%.3f, a2=%.3f, tau1=%.3f, tau2=%.3f', ...
        %     c, order2,idx, avglife2(idx), LT2Decay(idx), W2(1), W2(2), T2(1), T2(2)))
        %
        % % plot tri-exponential
        % nexttile
        % plot(WF(:,idx),'b.','MarkerSize',10)
        % hold on
        % plot(fit3(:,idx),'r-','LineWidth',1)
        % plot(res3(:,idx)-0.1,'m.')
        % yline(-0.1+0.025,'r--','res=0.025')
        % yline(-0.1-0.025,'r--')
        % text(150, -0.15, ['\Sigma res = ', num2str(res3SE)])
        % yticks((-0.2:0.1:1))
        % yticklabels({'','res=0','0','0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','1.0',})
        % grid on
        % grid minor
        % axis tight
        % ylim([-0.2 1])
        % legend('raw data', 'fitting', 'residue')
        % title(sprintf('Channel %d, Exponential fitting, order = %d, WF index %d,\n lifetime = %.3f(formula), lifetime = %.3f(decay),\n a1=%.3f, a2=%.3f, a3=%.3f, tau1=%.3f, tau2=%.3f, tau3=%.3f', ...
        %     c, order3,idx, avglife3(idx), LT3Decay(idx), W3(1), W3(2), W3(3), T3(1), T3(2), T3(3)))
        %
        % % plot quad-exponential
        % nexttile
        % plot(WF(:,idx),'b.','MarkerSize',10)
        % hold on
        % plot(fit4(:,idx),'r-','LineWidth',1)
        % plot(res4(:,idx)-0.1,'m.')
        % yline(-0.1+0.025,'r--','res=0.025')
        % yline(-0.1-0.025,'r--')
        % text(150, -0.15, ['\Sigma res = ', num2str(res4SE)])
        % yticks((-0.2:0.1:1))
        % yticklabels({'','res=0','0','0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','1.0',})
        % grid on
        % grid minor
        % axis tight
        % ylim([-0.2 1])
        % legend('raw data', 'fitting', 'residue')
        % title(sprintf('Channel %d, Exponential fitting, order = %d, WF index %d,\n lifetime = %.3f(formula), lifetime = %.3f(decay),\n a1=%.3f, a2=%.3f, a3=%.3f, a4=%.3f, tau1=%.3f, tau2=%.3f, tau3=%.3f, tau4=%.3f', ...
        %     c, order4, idx, avglife4(idx), LT4Decay(idx), W4(1), W4(2), W4(3), W4(4), T4(1), T4(2), T4(3), T4(4)))
        
    end
    
    exportgraphics(gcf,fullfile(DeConpath,[name,'Ch2.tiff']),'Resolution',600)
end
% close all
