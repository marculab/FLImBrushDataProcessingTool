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

for k = 1:length(DeConfile)
%% load in file
P = fullfile(DeConpath,DeConfile{k});
load(P)

figure('units','normalized','outerposition',[0 0 1 1])
tiledlayout(3,3)

for c = 1:3
%% set up channel
dataObj = Data{c}; %% set up your channel here
%% get Laguerre result
WFAligned = dataObj.channeldata.data;
SNR = dataObj.channeldata.SNR;
SNRIdx = SNR<30;
WFAligned(:,SNRIdx) = NaN;
WF = mean(WFAligned,2,'omitnan');
irf = dataObj.channeldata.iIRF;

% figure
% plot(WF); hold on;
% plot(irf)
% hold off


%% bi-exponential fit
order2 = 2;
[A2, T2, avglife2, intensity2, fit2, raw2, decay2] = multiexp_fit(WF,0.08,irf,order2,[]);
res2 = WF-fit2;
res2SE = sum(res2.^2);
LT2Decay = h_lifet(decay2,0.08);
W2 = A2./sum(A2);
% WTemp = WTemp';

%% tri-exponential fit
order3 = 3;
[A3, T3, avglife3, intensity3, fit3, raw3, decay3] = multiexp_fit(WF,0.08,irf,order3,[]);
res3 = WF-fit3;
res3SE = sum(res3.^2);
LT3Decay = h_lifet(decay3,0.08);
W3 = A3./sum(A3);

%% tri-exponential fit
order4 = 4;
[A4, T4, avglife4, intensity4, fit4, raw4, decay4] = multiexp_fit(WF,0.08,irf,order4,[]);
res4 = WF-fit4;
res4SE = sum(res4.^2);
LT4Decay = h_lifet(decay4,0.08);
W4 = A4./sum(A4);

%% plot Laguerre
idx = 1;

% plot bi-exponential
nexttile
plot(WF,'b.','MarkerSize',10)
hold on
plot(fit2,'r-','LineWidth',1)
plot(res2-0.1,'m.')
yline(-0.1+0.025,'r--','res=0.025')
yline(-0.1-0.025,'r--')
text(150, -0.15, ['\Sigma res = ', num2str(res2SE)])
yticks((-0.2:0.1:1))
yticklabels({'','res=0','0','0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','1.0',})
grid on
grid minor
axis tight
ylim([-0.2 1])
legend('raw data', 'fitting', 'residue')
title(sprintf('Channel %d, Exponential fitting, order = %d, WF index %d,\n lifetime = %.3f(formula), lifetime = %.3f(decay),\n a1=%.3f, a2=%.3f, tau1=%.3f, tau2=%.3f', ...
    c, order2,idx, avglife2(idx), LT2Decay(idx), W2(1), W2(2), T2(1), T2(2)))

% plot tri-exponential
nexttile
plot(WF,'b.','MarkerSize',10)
hold on
plot(fit3,'r-','LineWidth',1)
plot(res3-0.1,'m.')
yline(-0.1+0.025,'r--','res=0.025')
yline(-0.1-0.025,'r--')
text(150, -0.15, ['\Sigma res = ', num2str(res3SE)])
yticks((-0.2:0.1:1))
yticklabels({'','res=0','0','0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','1.0',})
grid on
grid minor
axis tight
ylim([-0.2 1])
legend('raw data', 'fitting', 'residue')
title(sprintf('Channel %d, Exponential fitting, order = %d, WF index %d,\n lifetime = %.3f(formula), lifetime = %.3f(decay),\n a1=%.3f, a2=%.3f, a3=%.3f, tau1=%.3f, tau2=%.3f, tau3=%.3f', ...
    c, order3,idx, avglife3(idx), LT3Decay(idx), W3(1), W3(2), W3(3), T3(1), T3(2), T3(3)))

% plot quad-exponential
nexttile
plot(WF,'b.','MarkerSize',10)
hold on
plot(fit4,'r-','LineWidth',1)
plot(res4-0.1,'m.')
yline(-0.1+0.025,'r--','res=0.025')
yline(-0.1-0.025,'r--')
text(150, -0.15, ['\Sigma res = ', num2str(res4SE)])
yticks((-0.2:0.1:1))
yticklabels({'','res=0','0','0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','1.0',})
grid on
grid minor
axis tight
ylim([-0.2 1])
legend('raw data', 'fitting', 'residue')
title(sprintf('Channel %d, Exponential fitting, order = %d, WF index %d,\n lifetime = %.3f(formula), lifetime = %.3f(decay),\n a1=%.3f, a2=%.3f, a3=%.3f, a4=%.3f, tau1=%.3f, tau2=%.3f, tau3=%.3f, tau4=%.3f', ...
    c, order4, idx, avglife4(idx), LT4Decay(idx), W4(1), W4(2), W4(3), W4(4), T4(1), T4(2), T4(3), T4(4)))

end
[filepath,name,ext]=fileparts(DeConfile{k});
exportgraphics(gcf,fullfile(DeConpath,[name,'run3.tiff']),'Resolution',600)

end
close all
