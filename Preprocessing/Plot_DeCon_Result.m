function handles = Plot_DeCon_Result(handles)
% get channel, line and point number
% update point slide if line slider is moved, otherwise do not update point
% slider
% clear all axis before plotting
axes(handles.axes_decay)
try
delete(handles.zoom_axes)
catch
end
    
cla reset;
hold off
axis on
box on
axes(handles.axes_residuals)
cla reset;
hold off
axis on
box on
axes(handles.axes_autocorr)
cla reset;
hold off
axis on
box on
axes(handles.axes_norm)
cla reset;
hold off
axis on
box on

Plot_option = getappdata(handles.figure_main,'Plot_option');
channel_choice = Plot_option.dataSource;
line_idx = Plot_option.lineIndex;
point_idx = Plot_option.pointIndex;
log_choice = Plot_option.logScale;
all_data = getappdata(handles.figure_main,'all_data');
raw_data_cell = getappdata(handles.figure_main,'raw_data_cell');
if (channel_choice>0 && channel_choice<5)
    if line_idx <= 0
        line_idx=1;
    end
    if line_idx>size(all_data,1)
        line_idx=size(all_data,1);
    end
    if point_idx<=0
        point_idx=1;
    end
    data_struct = all_data{line_idx,channel_choice};
    if isempty(data_struct)
        %     msgbox('No date for this line','modal');
        rawdata_to_plot = zeros(1,680);
        fitting_to_plot = zeros(1,680);
        residue_to_plot = zeros(1,680);
        autocorr_to_plot = zeros(1,680);
        decay_to_plot = zeros(1,680);
        iRF_to_plot = zeros(1,250);
        life_time = NaN;
        %     handles.fitting_to_plot = fitting;
        %     handles.data_struct_to_plot = data_struct;
        %     handles.residue = residue;
        axes(handles.axes_decay)
        title_temp = sprintf('No Data to display for Channel %2d, Line %2d.',channel_choice,line_idx);
        title(title_temp);
    else
        channeldata_all = get(data_struct,'channeldata');
        channeldata_all.data = channeldata_all.data';
        gain_list = channeldata_all.gain_list(:,1);
        decay_length = size(channeldata_all.data,2);
        dt = channeldata_all.dt;
        %     decon_idx = channeldata_all.decon_idx;
        % creat axis in ns for plot
        x_axis = dt*[1:decay_length];
        % iRF need a seperate x axis as it is usually only 250 point
        iRF_axis = dt*[1:length(channeldata_all.iIRF)];

        iRF_max = max(channeldata_all.iIRF);
        fitting = get(data_struct,'fit')';
        residue = get(data_struct,'res')';
        decay_to_plot = get(data_struct,'decay')';
        rawColor = [0 0 0];
        % temp code to scale the plotted data

        % check weather the point is deconvolved
        sudo_idx = find(ismember(channeldata_all.decon_idx,point_idx)); %deconvolved index, point_index in the real index.
        % if the point is deconvolved, get data
        if sudo_idx
            try
            gof = data_struct.stat_test.autoCorr.h(sudo_idx); % goodness of fit to filter data
            catch
                gof=1;
            end
            if gof
                rawColor = [0 0 1]; % if fit is good, use blue
            else
                rawColor = [0.5 0.5 0.5];% if fit is not good, use grey
            end
            data_max = max(channeldata_all.data(sudo_idx,:));
            rawdata_to_plot = channeldata_all.data(sudo_idx,:);
            factor = data_max/iRF_max;
            iRF_to_plot = channeldata_all.iIRF*factor;
            fitting_to_plot = fitting(sudo_idx,:);
            residue_to_plot = residue(sudo_idx,:);
            decay_to_plot = decay_to_plot(sudo_idx,:);
            life_time = data_struct.LTs(sudo_idx);
            gain = gain_list(sudo_idx);
            try
            autocorr_to_plot = data_struct.stat_test.autoCorr.autoCorrCurve{sudo_idx};
            catch
                autocorr_to_plot = zeros(1,680);
            end
            title_temp = sprintf('Channel %2d, Line %2d, Point %2d, Lifetime %3.2f, Gain %4.0f',channel_choice,line_idx,point_idx,life_time,gain);

            try
                figure(handles.h1);
            catch
                handles.h1 = figure('Position',[1 540 700 300]);
                movegui(handles.h1,'northeast');
            end

            raw_data = raw_data_cell{line_idx};
            BG = get(raw_data,'bg');
            raw_waveform = get(raw_data,'rawdata');
            bg_removed_data = get(raw_data,'data');
            %     figure(handles.h1)
            if ~isempty(raw_waveform)
            plot(raw_waveform(point_idx,:),'b','LineWidth', 2);
            end
            hold on
            if ~isempty(bg_removed_data)
            plot(bg_removed_data(sudo_idx,:),'r','LineWidth', 2);
            end
            if ~isempty(BG)
            plot(BG*raw_data.bg_scale_factor(sudo_idx),'g','LineWidth', 2);
            end
%           plot(raw_waveform(point_idx,:)-BG'*raw_data.bg_scale_factor(sudo_idx),'y','LineWidth', 2);
            hold off
            if ~isempty(BG)
                title(sprintf('Line %2d, Point %2d, Gain %4.0f, BG Scale factor %2.3f',line_idx,point_idx, gain, raw_data.bg_scale_factor(sudo_idx)));
            else
            title(sprintf('Line %2d, Point %2d, Gain %4.0f.',line_idx,point_idx, gain));
            end
            set(gca,'Position',[0.05 0.1 0.9 0.80])
            ylim([-0.05 0.5]);
            legend('Raw Data', 'BG Removed', 'Scaled BG','Location','northeast');
        else % if not set all plot to 0 excpet iRF
            data_max = 0.70;
            rawdata_to_plot = zeros(1,decay_length);
            factor = data_max/iRF_max;
            fitting_to_plot = zeros(1,decay_length);
            residue_to_plot = zeros(1,decay_length);
            autocorr_to_plot = zeros(1,680);
            decay_to_plot = zeros(1,decay_length);
            iRF_to_plot = channeldata_all.iIRF*factor;
            life_time = NaN;
            title_temp = sprintf('Channel %2d, Line %2d, Point %2d, data is either too low or saturated.',channel_choice,line_idx,point_idx);
        end
        %------------------------------------------------------------------
        % plot raw data and fitting----------------------------------------
        axes(handles.axes_decay)
        if log_choice
            semilogy(x_axis,rawdata_to_plot,'Color', rawColor,'LineWidth', 2)
            hold on
            semilogy(x_axis,fitting_to_plot,'r','LineWidth', 1);
            hold off
            ylim([0.000001 1])
            legend('Raw Data', 'Fitting');
        else
            plot(x_axis,rawdata_to_plot,'.-','Color', rawColor,'LineWidth', 1, 'MarkerSize',10);
            hold on
            plot(iRF_axis,iRF_to_plot,'g','LineWidth', 2);
            plot(x_axis,fitting_to_plot,'r','LineWidth', 1);
            hold off
            lgd = legend('Raw Data', 'iRF', 'Fitting');
            lgd.FontSize = 8;
        end
        ylabel('Voltage')
        title(title_temp);
        xlim([str2num(handles.x_min.String)  str2num(handles.x_max.String)]);
        ylim([-0.05 max(rawdata_to_plot)+0.05]);
        box on
        
        % plot of zoomed in of the rising edge-----------------------------
        handles.zoom_axes = axes('Position',[.55 .7 .2 .2]);
        box on
        plot(x_axis,rawdata_to_plot,'.-','Color', rawColor,'LineWidth', 1, 'MarkerSize',10);
        hold on
        plot(iRF_axis,iRF_to_plot,'g','LineWidth', 2);
        plot(x_axis,fitting_to_plot,'r','LineWidth', 1);
        hold off
        [~, max_idx] = max(iRF_to_plot);
        xlim([(0.08*max_idx-1) ceil(0.08*max_idx+1)])
        ylim([-0.02 max(rawdata_to_plot)*1.1])
        %--------------------------------------------------------------------------------------------------------------------
        % plot normalized residue. The result is in fraction
        axes(handles.axes_residuals)
        norm_residue_to_plot = residue_to_plot./data_max;
        plot(x_axis,norm_residue_to_plot,'r.-', 'LineWidth', 1,'MarkerSize',10);
        try
            ylim([-max(max(norm_residue_to_plot), abs(min(norm_residue_to_plot))) max(max(norm_residue_to_plot), abs(min(norm_residue_to_plot)))]);
        end
        %         ylim([-0.02 0.02]);
        ylabel('Residuals(a.u)');
        xlim([str2num(handles.x_min.String)  str2num(handles.x_max.String)])
        box on

        %---------------------------------------------------------------------------------------------------------------------
        % plot auto-correlation
        axes(handles.axes_autocorr)
        autocorr_to_plot = xcorr(residue_to_plot, 'coeff');
        plot(x_axis, autocorr_to_plot(length(x_axis):end), 'm', 'LineWidth', 2)
        % setup upper and lower bound of auto-correlation
        bounds(1) = 2 / sqrt(decay_length);
        bounds(2) = -bounds(1);
        hold on
        plot(x_axis, (bounds(1) * ones(size(x_axis))), 'k--')
        plot(x_axis, (bounds(2) * ones(size(x_axis))), 'k--')
        hold off
        ylabel('Res. Autocorrelation (a.u.)');
        try
            ylim([-max(max(autocorr_to_plot), abs(min(autocorr_to_plot))) max(max(autocorr_to_plot), abs(min(autocorr_to_plot)))]);
        end
        xlim([str2num(handles.x_min.String)  str2num(handles.x_max.String)])
        box on
        %------------------------------------------------------------------
        % plot fitted deccay-----------------------------------------------

        axes(handles.axes_norm)
        plot(x_axis, decay_to_plot/max(decay_to_plot), 'm', 'LineWidth', 2);
        ylabel('Norm. Decay (a.u.)');
        xlabel('ns')
        ylim([0 1])
        xlim([0 max(x_axis)])
        box on
        if sudo_idx
            tests = all_data{line_idx,channel_choice}.stat_test;
            if tests.chi2.h(sudo_idx)
                msg = sprintf('Chi2 TRUE - Data point passed Chi2 test. Chi2 value is %.2f.\n',tests.chi2.stat(sudo_idx));
                fprintf(msg);
            else
                msg = sprintf('Chi2 FALSE - Data point passed Chi2 test. Chi2 value is %.2f.\n',tests.chi2.stat(sudo_idx));
                fprintf(msg);
            end
            
            if tests.lbq.h(sudo_idx)
                display('LBQ TRUE - Data point passed Ljung Box Q test.');
            else
                display('LBQ FALSE - Data point did not Ljung Box Q test.');
            end
        else
            display('Not test result: data not deconvolved.');
        end
        
    end
end

drawnow
