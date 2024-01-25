function snrVSlt = plotsnrVSlt(Ch1DataObj, Ch2DataObj, Ch3DataObj, nbinss, snrlims)
    
    snrVSlt = figure;
    tiledlayout(3,1);

    % Plotting for Channel 1
    nexttile;
    histogram2(Ch1DataObj.Lg_LTs, Ch1DataObj.SNR, nbinss, 'DisplayStyle', 'tile');
    title('ch1');
    ylabel('SNR');
    colorbar;
    xlim([0 20]);
    ylim(snrlims);

    % Plotting for Channel 2
    nexttile;
    histogram2(Ch2DataObj.Lg_LTs, Ch2DataObj.SNR, nbinss, 'DisplayStyle', 'tile');
    title('ch2');
    ylabel('SNR');
    colorbar;
    xlim([0 20]);
    ylim(snrlims);

    % Plotting for Channel 3
    nexttile;
    histogram2(Ch3DataObj.Lg_LTs, Ch3DataObj.SNR, nbinss, 'DisplayStyle', 'tile');
    title('ch3');
    ylabel('SNR');
    colorbar;
    xlim([0 20]);
    ylim(snrlims);
    xlabel('Lifetime (ns)');
end