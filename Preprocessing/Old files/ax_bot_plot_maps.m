function handles = ax_bot_plot_maps(handles)

traget_to_plot = [];

channel_choice = get(handles.data_source,'Value');
plot_choice = get(handles.options_bottom,'Value');
LT_map = getappdata(handles.figure1,'LT_map_interped');
INT_map = getappdata(handles.figure1,'INT_map_interped');
SNR_map = getappdata(handles.figure1,'SNR_map_interped');
switch plot_choice
    case 1
        traget_to_plot = LT_map{channel_choice};
        axes(handles.ax_bottom);
        imagesc(traget_to_plot)
        colormap jet
        colorbar
        axis equal
        clear traget_to_plot
        title_text = sprintf('Lifetime of channel %d',channel_choice);
        title(title_text)
    case 2
        traget_to_plot = INT_map{channel_choice};
        axes(handles.ax_bottom);
        imagesc(traget_to_plot)
        colormap jet
        colorbar
        axis equal
        clear traget_to_plot
        title_text = sprintf('Intensity of channel %d',channel_choice);
        title(title_text)
    case 3
        traget_to_plot = LT_map{channel_choice};
        traget_to_plot = LT_map{channel_choice};
        msgbox('Function Not Yet Available','modal');
        clear traget_to_plot
        title_text = sprintf('Lifetime of channel %d',channel_choice);
        title(title_text)
        clear traget_to_plot
    case 4
        traget_to_plot = LT_map{channel_choice};
        msgbox('Function Not Yet Available','modal');
        clear traget_to_plot
    case 5
        traget_to_plot = SNR_map{channel_choice};
        axes(handles.ax_bottom);
        imagesc(traget_to_plot)
        colormap jet
        colorbar
        axis equal
        clear traget_to_plot
        
drawnow


end