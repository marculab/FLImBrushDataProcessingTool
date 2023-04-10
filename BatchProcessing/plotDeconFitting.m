function f = plotDeconFitting(plotObj, channel,plotIdx)
rawColor = [0 0 0 0.5]; % if fit is good, use blue
rawDataToPlot = get(plotObj,'wf_aligned',plotIdx);
rawDataMax = max(rawDataToPlot);
fitToPlotLG = get(plotObj,'fit',plotIdx);
residueToPlotLG = rawDataToPlot-fitToPlotLG;

gain = plotObj.gain(plotIdx);
CtrlV = plotObj.CtrlV(plotIdx);
LT_LG =  plotObj.Lg_LTs(plotIdx);
irfIdx = plotObj.irfIdx(plotIdx);

irf = plotObj.APDObj.irfTNorm(:,irfIdx);

t = plotObj.dtUp*(1:length(rawDataToPlot))';
%             xMaxEditField.Value = max(t); % set maximum x aixs value after loading data
% plot fitting
f = figure;
tiledlayout(3,1)
nexttile
plot(t,rawDataToPlot,'.-','Color', rawColor,'LineWidth', 1, 'MarkerSize',10)
hold('on')
plot((1:length(irf))*plotObj.dtUp,rawDataMax/max(irf)*irf,'Color',[0 0 1 0.5],'LineWidth', 1)
plot(t,fitToPlotLG,'m','LineWidth', 1.2)
hold('off')

lgd = legend('Raw Data', 'iRF', 'LG fit');
lgd.FontSize = 8;
title_temp = sprintf('Point %2d, Laguerre Lifetime %3.2f, Gain %4.0f.',plotIdx,LT_LG,gain);
title_temp = [channel ' ' title_temp];
title(title_temp)
ylim([-0.4 max(max(rawDataToPlot),0)+0.1])



% plot normalized residue
nexttile
normresidueToPlotLG = residueToPlotLG./rawDataMax*100;
plot(t,normresidueToPlotLG,'r.-', 'LineWidth', 1,'MarkerSize',10);
hold on
yline(2.5,'k--')
yline(-2.5,'k--')
hold off
lgd = legend('Laguerre', 'Location','northeast');
lgd.FontSize = 8;
try
    ylim([-max(max(norm_residue_to_plot), abs(min(norm_residue_to_plot))) max(max(norm_residue_to_plot), abs(min(norm_residue_to_plot)))]);
end
% xlim([xMinEditField.Value xMaxEditField.Value])
title('Residual (%)')

%plot auto correlation
nexttile
autoCorrToPlotLG = xcorr(residueToPlotLG, 'coeff');
plot( t, autoCorrToPlotLG(length(t):end), 'r.-', 'LineWidth', 2)
hold on
bounds(1) = 2 / sqrt(length(rawDataToPlot));
bounds(2) = -bounds(1);
plot( t, (bounds(1) * ones(size(t))), 'k--')
plot( t, (bounds(2) * ones(size(t))), 'k--')
hold('off')
ylabel('Res. Autocorrelation (a.u.)');
try
    ylim([-max(max(autocorr_to_plot), abs(min(autocorr_to_plot))) max(max(autocorr_to_plot), abs(min(autocorr_to_plot)))]);
end
% xlim([xMinEditField.Value xMaxEditField.Value])
lgd = legend('Laguerre', 'Location','northeast');
lgd.FontSize = 8;
end