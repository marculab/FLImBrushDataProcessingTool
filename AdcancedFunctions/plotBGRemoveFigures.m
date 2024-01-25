function fBGremove = plotBGRemoveFigures(Ch1DataObj, Ch2DataObj, Ch3DataObj, plotStep, bgLow, bgHigh, DigitizerNoise, ChWidthInPoints)
    fBGremove = figure;
    tiledlayout(3, 2);

    % Plotting for Channel 1
     nexttile;
    temp = sortData(Ch1DataObj, 'ascend');
    tempMin = min(temp(:));
    if isnan(tempMin)
        tempMin = 0;
    end
    plot(Ch1DataObj.averagedData(:, 1:plotStep:end));
    xline(bgLow, 'r--', 'BG', 'LineWidth', 2);
    xline(bgHigh, 'r--', 'BG', 'LineWidth', 2);
    yline(DigitizerNoise, 'm--', 'LineWidth', 2);
    %ylim([tempMin - 0.1 2]);
    ylim([max(tempMin - 0.1, 0) 2]);
    title(['Raw Data - Ch1']);

    % Plotting for DC&BG removed data
    nexttile;
    temp = sortData(Ch1DataObj, 'ascend');
    tempMin = min(temp(:));
    if isnan(tempMin)
        tempMin=0;
    end
    [~, maxIdx] = max(temp);
    maxIdx = median(maxIdx);
    plot(temp(:, 1:plotStep:end));
    xline(maxIdx - 200 + 1, 'c--', 'Truncation', 'LineWidth', 2);
    xline(maxIdx + ChWidthInPoints - 200, 'c--', 'Truncation', 'LineWidth', 2);
    xline(bgLow, 'r--', 'LineWidth', 2);
    xline(bgHigh, 'r--', 'LineWidth', 2);
    yline(DigitizerNoise, 'm--', 'LineWidth', 2);
    ylim([tempMin - 0.1 2]);
    title(['DC&BG Removed - Ch1']);

    % Plotting for Channel 2
     nexttile;
    temp = sortData(Ch2DataObj, 'ascend');
    tempMin = min(temp(:));
    if isnan(tempMin)
        tempMin = 0;
    end
    plot(Ch2DataObj.averagedData(:, 1:plotStep:end));
    xline(bgLow, 'r--', 'BG', 'LineWidth', 2);
    xline(bgHigh, 'r--', 'BG', 'LineWidth', 2);
    yline(DigitizerNoise, 'm--', 'LineWidth', 2);
    %ylim([tempMin - 0.1 2]);
    ylim([max(tempMin - 0.1, 0) 2]);
    title(['Raw Data - Ch2 ']);

    % Plotting for DC&BG removed data
    nexttile;
    temp = sortData(Ch2DataObj, 'ascend');
    tempMin = min(temp(:));
    if isnan(tempMin)
        tempMin=0;
    end
    [~, maxIdx] = max(temp);
    maxIdx = median(maxIdx);
    plot(temp(:, 1:plotStep:end));
    xline(maxIdx - 200 + 1, 'c--', 'Truncation', 'LineWidth', 2);
    xline(maxIdx + ChWidthInPoints - 200, 'c--', 'Truncation', 'LineWidth', 2);
    xline(bgLow, 'r--', 'LineWidth', 2);
    xline(bgHigh, 'r--', 'LineWidth', 2);
    yline(DigitizerNoise, 'm--', 'LineWidth', 2);
    ylim([tempMin - 0.1 2]);
    title(['DC&BG Removed - Ch2']);

    % Plotting for Channel 3
     nexttile;
    temp = sortData(Ch3DataObj, 'ascend');
    tempMin = min(temp(:));
    if isnan(tempMin)
        tempMin = 0;
    end
    plot(Ch3DataObj.averagedData(:, 1:plotStep:end));
    xline(bgLow, 'r--', 'BG', 'LineWidth', 2);
    xline(bgHigh, 'r--', 'BG', 'LineWidth', 2);
    yline(DigitizerNoise, 'm--', 'LineWidth', 2);
    %ylim([tempMin - 0.1 2]);
    ylim([max(tempMin - 0.1, 0) 2]);
    title(['Raw Data - Ch3']);

    % Plotting for DC&BG removed data
    nexttile;
    temp = sortData(Ch3DataObj, 'ascend');
    tempMin = min(temp(:));
    if isnan(tempMin)
        tempMin=0;
    end
    [~, maxIdx] = max(temp);
    maxIdx = median(maxIdx);
    plot(temp(:, 1:plotStep:end));
    xline(maxIdx - 200 + 1, 'c--', 'Truncation', 'LineWidth', 2);
    xline(maxIdx + ChWidthInPoints - 200, 'c--', 'Truncation', 'LineWidth', 2);
    xline(bgLow, 'r--', 'LineWidth', 2);
    xline(bgHigh, 'r--', 'LineWidth', 2);
    yline(DigitizerNoise, 'm--', 'LineWidth', 2);
    ylim([tempMin - 0.1 2]);
    title(['DC&BG Removed - Ch3']);

end

%%
function plotChannelFigures(ChannelDataObj, plotStep, bgLow, bgHigh, DigitizerNoise, ChWidthInPoints)
    % Plotting for raw data

    
    nexttile;
    temp = sortData(ChannelDataObj, 'ascend');
    tempMin = min(temp(:));
    if isnan(tempMin)
        tempMin = 0;
    end
    plot(ChannelDataObj.averagedData(:, 1:plotStep:end));
    xline(bgLow, 'r--', 'BG', 'LineWidth', 2);
    xline(bgHigh, 'r--', 'BG', 'LineWidth', 2);
    yline(DigitizerNoise, 'm--', 'LineWidth', 2);
    %ylim([tempMin - 0.1 2]);
    ylim([max(tempMin - 0.1, 0) 2]);
    title(['Raw Data - ']);

    % Plotting for DC&BG removed data
    nexttile;
    temp = sortData(ChannelDataObj, 'ascend');
    tempMin = min(temp(:));
    if isnan(tempMin)
        tempMin=0;
    end
    [~, maxIdx] = max(temp);
    maxIdx = median(maxIdx);
    plot(temp(:, 1:plotStep:end));
    xline(maxIdx - 200 + 1, 'c--', 'Truncation', 'LineWidth', 2);
    xline(maxIdx + ChWidthInPoints - 200, 'c--', 'Truncation', 'LineWidth', 2);
    xline(bgLow, 'r--', 'LineWidth', 2);
    xline(bgHigh, 'r--', 'LineWidth', 2);
    yline(DigitizerNoise, 'm--', 'LineWidth', 2);
    ylim([tempMin - 0.1 2]);
    title(['DC&BG Removed - ']);
end