function estimatebkgd(Ch1DataObj, Ch2DataObj, Ch3DataObj, runName, saveFigFullPath)
    % Create a new figure window and tiled layout
    fCh1 = figure;                             
    tiledlayout('flow');                      

    % Process each channel
    processChannel(Ch1DataObj, 'Channel 1');
    processChannel(Ch2DataObj, 'Channel 2');
    processChannel(Ch3DataObj, 'Channel 3');

    % Save figure as .mat file
    saveFileName_fig = strcat(runName, '_Fig_BkgdEstimation.fig');
    saveFigFullPath = fullfile(saveFigFullPath, saveFileName_fig);
    savefig(fCh1, saveFigFullPath, 'compact');

    % Close the figure to clean up
    close(fCh1);


    % Define a helper function to process each channel
    function processChannel(DataObj, channelName)
        bg = DataObj.bg;
        bgLow = DataObj.bgLow;
        bgHigh = DataObj.bgHigh;

        % Plot the raw data from DataObj
        nexttile
        plot(DataObj.rawDataUpsampled(:,1:500:end)) % selecting every 500th column of the raw data matrix
        hold on;
        plot(bg, 'r', 'LineWidth', 2);
        xline(bgLow, 'b--', 'LineWidth', 1.5); 
        xline(bgHigh, 'b--', 'LineWidth', 1.5); 
        xlim([500 2000]);
        title([channelName ' Upsampled Waveforms and measured Bkgd Waveform']);

        nexttile
        plot(DataObj.preProcessedData(:,1:500:end));
        title([channelName ' BG removed Data'])

        tempMax = max(DataObj.rawDataUpsampled);
        tempIdx = find(tempMax > 0.9 & tempMax < 1.1);

        if isempty(tempIdx)
            [~, I] = sort(tempMax, 'descend');
            estimatedBG = DataObj.rawDataUpsampled(:, I(1:100));
        else
            temp = DataObj.rawDataUpsampled(:, tempIdx);
            tempCtrlV = DataObj.rawCtrlV(tempIdx);
            [~, I] = sort(tempCtrlV, 'descend');
            estimatedBG = temp(:, I(1:100));
        end
        estimatedBGAvg = nanmean(estimatedBG, 2);


        nexttile
        plot(estimatedBG)
        hold on
        xline(bgLow, 'b--', 'LineWidth', 1.5);
        xline(bgHigh, 'b--', 'LineWidth', 1.5); 
        h = plot(DataObj.bg, 'r-', 'LineWidth', 2);
        h2 = plot(estimatedBGAvg, 'g-', 'LineWidth', 2);
        ylim([-0.2 1.4])
        title([channelName ' Estimated BG from raw data'])
        legend([h, h2], 'Measured BG', 'Estimated BG')
        title([channelName ' Upsampled Waveforms and Bkgd Waveforms'])

        % Calculate average signal level for fiber background range for each
        % waveform and histogram
        avgSignalLevel = nanmean(DataObj.rawDataUpsampled(680:980, :), 1);

        mean_bg_est = nanmean(estimatedBGAvg(680:980, :), 1);
        mean_bg_meas = nanmean(bg(680:980, :), 1);
        per_5 = prctile(avgSignalLevel, 95);
        per_10 = prctile(avgSignalLevel, 90);

        % Visualization - Compare distributions
        nexttile
        histogram(avgSignalLevel, 100, 'Normalization', 'probability');
        hold on;
        xline(per_5, 'r', 'LineWidth', 1.5);
        xline(per_10, 'r', 'LineWidth', 1.5);
        h = xline(mean_bg_est, 'g', 'LineWidth', 1.5);
        h2 = xline(mean_bg_meas, 'b', 'LineWidth', 1.5);
        hold off;
        title([channelName ' Distribution of Average Signal Level (100-200)']);
        xlabel('Average Signal Level');
        ylabel('Probability');
        legend([h, h2], 'Mean est bg', 'Mean meas bg');

        % Check condition to replace background
        if mean_bg_meas < per_5
            DataObj.bg = estimatedBGAvg; % overwrite existing BG with estimated BG
            DataObj.BgEstimated = 1;     % set BG estimation flag
        end
    end
end