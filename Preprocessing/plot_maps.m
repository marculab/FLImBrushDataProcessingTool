function handles = plot_maps(handles,axes_choice)
Image_option = getappdata(handles.figure_main,'Image_option');
Plot_option = getappdata(handles.figure_main,'Plot_option');

if strcmp(axes_choice,'top')
    axes(handles.axes_top);
    plot_choice = Image_option.imageType(1);
else
    axes(handles.axes_bottom);
    plot_choice = Image_option.imageType(2);
end
target_to_plot = [];

S_Scan_flag = Image_option.Sscan;
pixels_to_shift = Image_option.pixelShift;

X_unit = Image_option.aspectRatio(1);
Y_unit = Image_option.aspectRatio(2);
Z_unit = 1;

INT_sum = 0;

channel_choice = Plot_option.dataSource;

if handles.Interp_check_box.Value %check wether to plot interploted map
    LT_map = getappdata(handles.figure_main,'LT_map_interped');
    INT_map = getappdata(handles.figure_main,'INT_map_interped');
    INT_R_map = getappdata(handles.figure_main,'INT_R_map_interped');
    SNR_map = getappdata(handles.figure_main,'SNR_map_interped');
    for i=1:length(INT_map)
        if ~isempty(INT_map{i})
    INT_sum =INT_sum+INT_map{i};
        end
    end
else
    LT_map = getappdata(handles.figure_main,'LT_map');
    INT_map = getappdata(handles.figure_main,'INT_map');
    INT_R_map = getappdata(handles.figure_main,'INT_R_map');
    SNR_map = getappdata(handles.figure_main,'SNR_map');
    for i=1:length(INT_map)
        if ~isempty(INT_map{i})
    INT_sum =INT_sum+INT_map{i};
        end
    end
end

SNR_threshold = str2double(handles.SNR_Threshold.String);
SNR_mask = ones(size(SNR_map{channel_choice}));
SNR_mask(SNR_map{channel_choice}<SNR_threshold)=0;

lt_low = str2double(handles.LT_low.String);
lt_high = str2double(handles.LT_high.String);



switch plot_choice
    case 1
        target_to_plot = LT_map{channel_choice}.*SNR_mask;
        if S_Scan_flag
            target_to_plot = flip_even_rows(target_to_plot,pixels_to_shift);
        end
        target_to_plot = flipud(target_to_plot);
        imagesc(target_to_plot,[lt_low lt_high]);
        daspect([X_unit Y_unit Z_unit])
        colormap jet
        colorbar
        clear target_to_plot
        
        title_text = sprintf('Lifetime of channel %d',channel_choice);
        title(title_text)
    case 2
        target_to_plot = INT_map{channel_choice}.*SNR_mask;
        if S_Scan_flag
            target_to_plot = flip_even_rows(target_to_plot,pixels_to_shift);
        end
        target_to_plot = flipud(target_to_plot);
        imagesc(target_to_plot)
        daspect([X_unit Y_unit Z_unit])
        colormap jet
        colorbar
        clear target_to_plot
        title_text = sprintf('Abs Intensity of channel %d',channel_choice);
        title(title_text)
    case 3
        target_to_plot = LT_map{channel_choice};
        target_to_plot = flipud(target_to_plot);
        msgbox('Function Not Yet Available','modal');
        clear target_to_plot
        title_text = sprintf('Lifetime of channel %d',channel_choice);
        title(title_text)
    case 4
        target_to_plot = INT_R_map{channel_choice}.*SNR_mask;
        if S_Scan_flag
            target_to_plot = flip_even_rows(target_to_plot,pixels_to_shift);
        end
        target_to_plot = flipud(target_to_plot);
        imagesc(target_to_plot)
        daspect([X_unit Y_unit Z_unit])
        colormap jet
        colorbar
        clear target_to_plot
        title_text = sprintf('AUC of channel %d',channel_choice);
        title(title_text)
    case 5
        target_to_plot = SNR_map{channel_choice};
        if S_Scan_flag
            target_to_plot = flip_even_rows(target_to_plot,pixels_to_shift);
        end
        target_to_plot = flipud(target_to_plot);
        imagesc(target_to_plot)
        daspect([X_unit Y_unit Z_unit])
        colormap jet
        colorbar
        
        %         axis equal
        clear target_to_plot
        title_text = sprintf('SNR of channel %d',channel_choice);
        title(title_text)
        
    case 6
        traget_to_plot = Taus_1_map{channel_choice}.*SNR_mask;
        if S_Scan_flag
            traget_to_plot = flip_even_rows(traget_to_plot,pixels_to_shift);
        end
        target_to_plot = flipud(target_to_plot);
        imagesc(traget_to_plot,[0 5])
        daspect([X_unit Y_unit Z_unit])
        colormap jet
        colorbar
        %         axis equal
        clear traget_to_plot
        title_text = sprintf('SNR of channel %d',channel_choice);
        title(title_text)
    case 7
        traget_to_plot = Taus_2_map{channel_choice}.*SNR_mask;
        if S_Scan_flag
            traget_to_plot = flip_even_rows(traget_to_plot,pixels_to_shift);
        end
        target_to_plot = flipud(target_to_plot);
        imagesc(traget_to_plot,[5 20])
        daspect([X_unit Y_unit Z_unit])
        colormap jet
        colorbar
        %         axis equal
        clear traget_to_plot
        title_text = sprintf('SNR of channel %d',channel_choice);
        title(title_text)
    case 8
        traget_to_plot = Taus_3_map{channel_choice}.*SNR_mask;
        
        if S_Scan_flag
            target_to_plot = flip_even_rows(target_to_plot,pixels_to_shift);
        end
        target_to_plot = flipud(target_to_plot);
        imagesc(target_to_plot)
        daspect([X_unit Y_unit Z_unit])
        colormap jet
        colorbar
        
        %         axis equal
        clear traget_to_plot
        title_text = sprintf('SNR of channel %d',channel_choice);
        title(title_text)
        
    case 9
        traget_to_plot = Weight_1_map{channel_choice}.*SNR_mask;
        if S_Scan_flag
            traget_to_plot = flip_even_rows(traget_to_plot,pixels_to_shift);
        end
        target_to_plot = flipud(target_to_plot);
        imagesc(traget_to_plot,[0 1])
        daspect([X_unit Y_unit Z_unit])
        colormap jet
        colorbar
        %         axis equal
        clear traget_to_plot
        title_text = sprintf('SNR of channel %d',channel_choice);
        title(title_text)
        
    case 10
        traget_to_plot = Weight_2_map{channel_choice}.*SNR_mask;
        if S_Scan_flag
            traget_to_plot = flip_even_rows(traget_to_plot,pixels_to_shift);
        end
        target_to_plot = flipud(target_to_plot);
        imagesc(traget_to_plot, [0 1])
        daspect([X_unit Y_unit Z_unit])
        colormap jet
        colorbar
        %         axis equal
        clear traget_to_plot
        title_text = sprintf('SNR of channel %d',channel_choice);
        title(title_text)
    case 11
        traget_to_plot = Weight_3_map{channel_choice}.*SNR_mask;
        if S_Scan_flag
            traget_to_plot = flip_even_rows(traget_to_plot,pixels_to_shift);
        end
        target_to_plot = flipud(target_to_plot);
        imagesc(traget_to_plot, [0 1])
        daspect([X_unit Y_unit Z_unit])
        colormap jet
        colorbar
        %         axis equal
        clear traget_to_plot
        title_text = sprintf('SNR of channel %d',channel_choice);
        title(title_text)
        
end
cmap = jet(256);
cmap(1,:) = zeros(1,3);
colormap(cmap)
drawnow

    function map_out = flip_even_rows(map_in,shift)
        map_out = map_in;
        map_out(2:2:end,:) = circshift(fliplr(map_out(2:2:end,:)),shift,2);
    end

end
