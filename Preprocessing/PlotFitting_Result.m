function PlotFitting_Result(handles,point_idx)

axes(handles.ax_decay)
channeldata_all = get(handles.laguerre_model_to_plot,'channeldata');
channeldata_all.data = channeldata_all.data';
data_max = max(channeldata_all.data(point_idx,:));
iRF_max = max(channeldata_all.iIRF);
factor = data_max/iRF_max;
plot(channeldata_all.data(point_idx,:),'b','LineWidth', 1)
hold on
plot(channeldata_all.iIRF.*factor,'g');
plot(handles.fitting_to_plot(point_idx,:),'r');
hold off
xlabel('Point')
ylabel('Voltage')

axes(handles.ax_residuals)
plot(handles.residue(point_idx,:),'r', 'LineWidth', 3);
ylim([-5e-3 5e-3]);
ylabel('Residuals(a.u)');

end

