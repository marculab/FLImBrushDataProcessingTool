% Define the function
function plotAndSaveFig(DataObj, bgLow, bgHigh, runName, saveFigFullPath, ChX)
    % Create a new figure window and tiled layout
    fCh1 = figure;                             % Create a new figure window
    tiledlayout('flow');                       % Create a tiled layout for subplots
        
    % Plot the raw data from DataObj
    nexttile
    plot(DataObj.rawData(:,:)) % selecting every 500th column of the raw data matrix
    grid on
    title('Raw Data')
    
    % % Plot the raw data with DC offset removed
    % nexttile
    % plot(Ch1DataObj.rawDataDCRemoved(:,1:500:end)); % selecting every 500th column of the raw data matrix 
    % title('DC removed Raw Data')
    
    % Plot the upsampled data
    % nexttile
    % plot(Ch1DataObj.rawDataUpsampled(:,1:500:end)); % selecting every 500th column of the raw data matrix 
    % title('Upsampled Data')
    
    
    nexttile
    plot(DataObj.averagedData(:,:)); hold on;
    plot(DataObj.bg, 'r', 'LineWidth', 2);
    xline(bgLow, 'b--', 'LineWidth', 1.5); 
    xline(bgHigh, 'b--', 'LineWidth', 1.5); 
    xlim([500 2000]);
    grid on
    title('Averaged Waveforms and Background Waveform');
    
    %nexttile
    %plot(Ch1DataObj.preProcessedData(:,1:500:end));
    %title('BG removed Data')
    

    nexttile
    plot(DataObj.dataT(:,:));
    grid on
    title('Truncated Data')
    

    % Plot lifetime values with gain > 300 in red and gain < 300 in blue
    Gain_Thr = 300;
    LT = DataObj.Lg_LTs;

    nexttile(4)
    hold on;
    scatter(find(DataObj.gain > Gain_Thr), LT(DataObj.gain > Gain_Thr), 10, 'r.');
    scatter(find(DataObj.gain <= Gain_Thr), LT(DataObj.gain <= Gain_Thr), 10, 'b.');
    title('Lifetime Values by Gain');
    legend('Gain > 300', 'Gain <= 300');
    xlabel('Index');
    ylabel('Lifetime');
    hold off;
    
    
    % Plot gain values
    nexttile
    plot(DataObj.gain);
    title('Gain Values');
    xlabel('Index');
    ylabel('Gain');
    
        
    % Plot lifetime values with SNR < 40 in red and SNR > 40 in blue
    % SNR_Thr = 40;
    % nexttile
    % hold on;
    % scatter(find(DataObj.SNR < SNR_Thr), DataObj.Lg_LTs(DataObj.SNR < SNR_Thr), 20,'r.');
    % scatter(find(DataObj.SNR >= SNR_Thr),DataObj.Lg_LTs(DataObj.SNR >= SNR_Thr), 20, 'b.');
    % title('Lifetime Values by SNR');
    % legend('SNR < 40', 'SNR >= 40');
    % xlabel('Index');
    % ylabel('Lifetime');
    % hold off;
    
    % % Plot SNR values
    nexttile
    plot(DataObj.SNR);
    title('SNR Values');
    xlabel('Index');
    ylabel('SNR');
    
    % 2D histogram of SNR vs. Gain
    nexttile
    histogram2(DataObj.Lg_LTs, DataObj.SNR, [500 500], 'DisplayStyle', 'tile');
    colorbar;
    title('2D Histogram of LT vs. SNR');
    xlabel('LT');
    ylabel('SNR');
    
    % Histogram of lifetime values
    nexttile
    histogram(DataObj.Lg_LTs);
    title('Histogram of Lifetime Values');
    xlabel('Lifetime');
    ylabel('Frequency');
    
    % % Histogram of intensity ratios
    % nexttile
    % histogram(Ch1DataObj.Lg_INTs);
    % title('Histogram of Intensity Ratios');
    % xlabel('Intensity Ratio');
    % ylabel('Frequency');
    
    % Phasor plot
    nexttile
    scatter(DataObj.Ph_H1S, DataObj.Ph_H1G, '.');
    title('Phasor Plot');
    xlabel('Phasor X');
    ylabel('Phasor Y');
    axis equal; % Ensure the plot is circular
    hold on;
    % Draw the half unit circle
    theta = linspace(0, pi, 100);
    x = cos(theta);
    y = sin(theta);
    plot(x, y, 'r--');
    hold off;
    
    % Save figure as .mat file
    saveFileName_fig = strcat(runName, '_', ChX, '_Figure.fig');
    [parentDir, ~, ~] = fileparts(saveFigFullPath);
    saveFigPath = fullfile(parentDir,saveFileName_fig);
    savefig(fCh1, saveFigPath, 'compact');

    % Save figures as jpg file
    set(fCh1, 'PaperUnits', 'millimeters');
    set(fCh1, 'PaperSize', [297 210]);
    set(fCh1, 'PaperPosition', [0 0 297 210]);
    print(fCh1, saveFigPath, '-djpeg', '-r300'); % 300 dpi resolution


    % close the figure to clean up
    close(fCh1);
end