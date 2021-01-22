function handles = load_processed_data(handles)

[name_cell,PathName] = uigetfile('*Decon.mat','Select the processed mat file','Multiselect', 'on',handles.rootFolder);


if ischar(name_cell)
    name_cell = {name_cell};
    Num_of_files = 1;
else
    name_cell  = name_cell';
    Num_of_files = length(name_cell);
end
% reinitilize everything to clear data
INT1_map  = cell(size(name_cell));
LT1_map  = cell(size(name_cell));
INT2_map  = cell(size(name_cell));
LT2_map  = cell(size(name_cell));
INT3_map  = cell(size(name_cell));
LT3_map  = cell(size(name_cell));
INT4_map  = cell(size(name_cell));
LT4_map  = cell(size(name_cell));

INT1_map_interped  = cell(size(name_cell));
LT1_map_interped  = cell(size(name_cell));
INT2_map_interped  = cell(size(name_cell));
LT2_map_interped  = cell(size(name_cell));
INT3_map_interped  = cell(size(name_cell));
LT3_map_interped  = cell(size(name_cell));
INT4_map_interped  = cell(size(name_cell));
LT4_map_interped  = cell(size(name_cell));

SNR1_map  = cell(size(name_cell));
SNR2_map  = cell(size(name_cell));
SNR3_map  = cell(size(name_cell));
SNR4_map  = cell(size(name_cell));

SNR1_map_interped  = cell(size(name_cell));
SNR2_map_interped  = cell(size(name_cell));
SNR3_map_interped  = cell(size(name_cell));
SNR4_map_interped  = cell(size(name_cell));
% main GUI variable to save all the data structure
all_data = cell(size(name_cell,1),4);

raw_data_cell = cell(size(name_cell,1),1);


for i=1:Num_of_files

    if  Num_of_files==1
        load(cell2mat(fullfile(PathName,name_cell)));
    else
    load(cell2mat(fullfile(PathName,name_cell(i))));
    end

    if  i == 1; setappdata(handles.figure_main,'n_pix',full_data_size(1));end %set number of pixels for 1st file loaded
    all_data(i,:) = Data;
    try
    raw_data_cell{i} = rawDataClass;
    catch
        
    end
    INT1_map{i}  = spec_int{1}';
    LT1_map{i}  = lifet_avg{1}';
    SNR1_map{i}  = SNR{1}';
    INT2_map{i} = spec_int{2}';
    LT2_map{i}  = lifet_avg{2}';
    SNR2_map{i}  = SNR{2}';
    INT3_map{i}  = spec_int{3}';
    LT3_map{i}  = lifet_avg{3}';
    SNR3_map{i}  = SNR{3}';
    INT4_map{i} = spec_int{4}';
    LT4_map{i}  = lifet_avg{4}';
    SNR4_map{i}  = SNR{4}';

    INT1_map_interped{i}  = spec_int_interped{1}';
    LT1_map_interped{i}  = lifet_avg_interped{1}';
    SNR1_map_interped{i}  = SNR_interped{1}';
    INT2_map_interped{i} = spec_int_interped{2}';
    LT2_map_interped{i}  = lifet_avg_interped{2}';
    SNR2_map_interped{i}  = SNR_interped{2}';
    INT3_map_interped{i}  = spec_int_interped{3}';
    LT3_map_interped{i}  = lifet_avg_interped{3}';
    SNR3_map_interped{i}  = SNR_interped{3}';
    INT4_map_interped{i} = spec_int_interped{4}';
    LT4_map_interped{i}  = lifet_avg_interped{4}';
    SNR4_map_interped{i}  = SNR_interped{4}';



end

INT_map{1} = cell2mat(INT1_map);
INT_map{2} = cell2mat(INT2_map);
INT_map{3} = cell2mat(INT3_map);
INT_map{4} = cell2mat(INT4_map);
INT_map_interped{1} = cell2mat(INT1_map_interped);
INT_map_interped{2} = cell2mat(INT2_map_interped);
INT_map_interped{3} = cell2mat(INT3_map_interped);
INT_map_interped{4} = cell2mat(INT4_map_interped);
setappdata(handles.figure_main,'INT_map',INT_map);
setappdata(handles.figure_main,'INT_map_interped',INT_map_interped);


LT_map{1} = cell2mat(LT1_map);
LT_map{2} = cell2mat(LT2_map);
LT_map{3} = cell2mat(LT3_map);
LT_map{4} = cell2mat(LT4_map);
LT_map_interped{1} = cell2mat(LT1_map_interped);
LT_map_interped{2} = cell2mat(LT2_map_interped);
LT_map_interped{3} = cell2mat(LT3_map_interped);
LT_map_interped{4} = cell2mat(LT4_map_interped);
setappdata(handles.figure_main,'LT_map',LT_map);
setappdata(handles.figure_main,'LT_map_interped',LT_map_interped);

SNR_map{1} = cell2mat(SNR1_map);
SNR_map{2} = cell2mat(SNR2_map);
SNR_map{3} = cell2mat(SNR3_map);
SNR_map{4} = cell2mat(SNR4_map);
SNR_map_interped{1} = cell2mat(SNR1_map_interped);
SNR_map_interped{2} = cell2mat(SNR2_map_interped);
SNR_map_interped{3} = cell2mat(SNR3_map_interped);
SNR_map_interped{4} = cell2mat(SNR4_map_interped);
setappdata(handles.figure_main,'SNR_map',SNR_map);
setappdata(handles.figure_main,'SNR_map_interped',SNR_map_interped);

setappdata(handles.figure_main,'all_data',all_data);
setappdata(handles.figure_main,'raw_data_cell',raw_data_cell);
