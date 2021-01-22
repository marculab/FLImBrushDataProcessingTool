function handles = pick_channel_data(handles,line_idx,channel_choice)
if isempty()
laguerre_model = handles.all_data{line_idx,channel_choice};
if ~isempty()
fitting = get(laguerre_model,'fit')';
residue = get(laguerre_model,'res')';
handles.fitting_to_plot = fitting;
handles.laguerre_model_to_plot = laguerre_model;
handles.residue = residue;
end

end
