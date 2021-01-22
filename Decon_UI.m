function varargout = Decon_UI(varargin)
% DECON_UI MATLAB code for Decon_UI.fig
%      DECON_UI, by itself, creates a new DECON_UI or raises the existing
%      singleton*.
%
%      H = DECON_UI returns the handle to a new DECON_UI or the handle to
%      the existing singleton*.
%
%      DECON_UI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DECON_UI.M with the given input arguments.
%
%      DECON_UI('Property','Value',...) creates a new DECON_UI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Decon_UI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Decon_UI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Decon_UI

% Last Modified by GUIDE v2.5 28-Jul-2017 19:00:30

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Decon_UI_OpeningFcn, ...
                   'gui_OutputFcn',  @Decon_UI_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before Decon_UI is made visible.
function Decon_UI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Decon_UI (see VARARGIN)


% save the handle of the main UI as 'figure_main' in current handles
movegui(handles.figure_decon,'west');
figure_main = varargin{1};
handles.figure_main = figure_main;
handles.View_decon = varargin{2};
% update Decon_UI controls according to main UI defaults
Files = getappdata(figure_main,'Files');
Decon_param = getappdata(figure_main,'Decon_param');
iIRF = getappdata(figure_main,'irf');

set(handles.uitable_files,'Data',[Files.selector Files.namelist]);

set(handles.edit_time_res,'String',sprintf('%g',Decon_param.time_res));
set(handles.edit_interp_res,'String',sprintf('%g',Decon_param.interpolation));
set(handles.edit_bw,'String',sprintf('%g',Decon_param.bw));
set(handles.edit_ch_width,'String',sprintf('%g',Decon_param.ch_width));
if isnan(Decon_param.LG_param.alpha)
    set(handles.edit_alpha,'String','auto');
else
    set(handles.edit_alpha,'String',sprintf('%g',Decon_param.LG_param.alpha));
end
set(handles.edit_LG_order,'String',sprintf('%g',Decon_param.LG_param.order));
set(handles.popupmenu_num_exp,'Value',Decon_param.Exp_param.model);
set(handles.edit_daq_time,'String',sprintf('%g',Decon_param.DAQ_time));
set(handles.edit_amp_low,'String',sprintf('%g',Decon_param.amplitude_window(1)));
set(handles.edit_amp_high,'String',sprintf('%g',Decon_param.amplitude_window(2)));
set(handles.checkbox_ch1,'Value',Decon_param.Channels(1));
set(handles.checkbox_ch2,'Value',Decon_param.Channels(2));
set(handles.checkbox_ch3,'Value',Decon_param.Channels(3));
set(handles.checkbox_ch4,'Value',Decon_param.Channels(4));
set(handles.popupmenu_decon_method,'Value',Decon_param.Decon_option);

set(handles.edit_iIRF_filename,'String',iIRF.name);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Decon_UI wait for user response (see UIRESUME)
% uiwait(handles.figure_decon);


% --- Outputs from this function are returned to the command line.
function varargout = Decon_UI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles;

% --- Executes during object deletion, before destroying properties.
function figure_decon_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to figure_decon (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.View_decon,'Checked','Off');
setappdata(handles.figure_main,'handles_decon',{});


function edit_iIRF_filename_Callback(hObject, eventdata, handles)
% hObject    handle to edit_iIRF_filename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_iIRF_filename as text
%        str2double(get(hObject,'String')) returns contents of edit_iIRF_filename as a double


% --- Executes during object creation, after setting all properties.
function edit_iIRF_filename_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_iIRF_filename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in button_select_all.
function button_select_all_Callback(hObject, eventdata, handles)
% hObject    handle to button_select_all (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Files = getappdata(handles.figure_main,'Files'); % check whether there is loaded data file
if ~isempty(Files)
    Files.selector = num2cell(logical(ones(size(Files.namelist))));
    set(handles.uitable_files,'Data',[Files.selector Files.namelist]);
end
guidata(hObject, handles);

% --- Executes on button press in button_remove_all.
function button_remove_all_Callback(hObject, eventdata, handles)
% hObject    handle to button_remove_all (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Files = struct('folderlist',{},'namelist',{},'selector',{});
setappdata(handles.figure_main,'Files',Files);
set(handles.uitable_files,'Data',[]);

guidata(hObject, handles);

% --- Executes on button press in button_remove_unselected.
function button_remove_unselected_Callback(hObject, eventdata, handles)
% hObject    handle to button_remove_unselected (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Files = getappdata(handles.figure_main,'Files');
if ~isempty(Files)
    selector = cell2mat(Files.selector(:,1));
    Files.namelist = Files.namelist(selector);
    Files.folderlist = Files.folderlist(selector);
    Files.selector = num2cell(logical(ones(size(Files.namelist))));
    setappdata(handles.figure_main,'Files',Files);
    set(handles.uitable_files,'Data',[Files.selector Files.namelist]);
    
end
guidata(hObject, handles);


% --- Executes on button press in button_unselect_all.
function button_unselect_all_Callback(hObject, eventdata, handles)
% hObject    handle to button_unselect_all (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Files = getappdata(handles.figure_main,'Files');
if ~isempty(Files)
    Files.selector = num2cell(logical(zeros(size(Files.namelist))));
    set(handles.uitable_files,'Data',[Files.selector Files.namelist]);
end
guidata(hObject, handles);


function edit_time_res_Callback(hObject, eventdata, handles)
% hObject    handle to edit_time_res (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_time_res as text
%        str2double(get(hObject,'String')) returns contents of edit_time_res as a double
Decon_param = getappdata(handles.figure_main,'Decon_param');
time_res_string = get(handles.edit_time_res,'String');
Decon_param.time_res = str2double(time_res_string);
setappdata(handles.figure_main,'Decon_param',Decon_param);

% --- Executes during object creation, after setting all properties.
function edit_time_res_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_time_res (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_interp_res_Callback(hObject, eventdata, handles)
% hObject    handle to edit_interp_res (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_interp_res as text
%        str2double(get(hObject,'String')) returns contents of edit_interp_res as a double
Decon_param = getappdata(handles.figure_main,'Decon_param');
interp_res_string = get(handles.edit_interp_res,'String');
Decon_param.interpolation = str2double(interp_res_string);
setappdata(handles.figure_main,'Decon_param',Decon_param);

% --- Executes during object creation, after setting all properties.
function edit_interp_res_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_interp_res (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_bw_Callback(hObject, eventdata, handles)
% hObject    handle to edit_bw (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_bw as text
%        str2double(get(hObject,'String')) returns contents of edit_bw as a double
Decon_param = getappdata(handles.figure_main,'Decon_param');
bw_string = get(handles.edit_bw,'String');
Decon_param.bw = str2double(bw_string);
setappdata(handles.figure_main,'Decon_param',Decon_param);

% --- Executes during object creation, after setting all properties.
function edit_bw_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_bw (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_ch_width_Callback(hObject, eventdata, handles)
% hObject    handle to edit_ch_width (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_ch_width as text
%        str2double(get(hObject,'String')) returns contents of edit_ch_width as a double
Decon_param = getappdata(handles.figure_main,'Decon_param');
ch_width_string = get(handles.edit_ch_width,'String');
Decon_param.ch_width = str2double(ch_width_string);
setappdata(handles.figure_main,'Decon_param',Decon_param);


% --- Executes during object creation, after setting all properties.
function edit_ch_width_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_ch_width (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_alpha_Callback(hObject, eventdata, handles)
% hObject    handle to edit_alpha (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_alpha as text
%        str2double(get(hObject,'String')) returns contents of edit_alpha as a double
Decon_param = getappdata(handles.figure_main,'Decon_param');
alpha_string = get(handles.edit_alpha,'String');
Decon_param.LG_param.alpha = str2double(alpha_string);
setappdata(handles.figure_main,'Decon_param',Decon_param);

% --- Executes during object creation, after setting all properties.
function edit_alpha_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_alpha (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_amp_low_Callback(hObject, eventdata, handles)
% hObject    handle to edit_amp_low (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_amp_low as text
%        str2double(get(hObject,'String')) returns contents of edit_amp_low as a double
Decon_param = getappdata(handles.figure_main,'Decon_param');
amp_low_string = get(handles.edit_amp_low,'String');
Decon_param.amplitude_window(1) = str2double(amp_low_string);
setappdata(handles.figure_main,'Decon_param',Decon_param);

% --- Executes during object creation, after setting all properties.
function edit_amp_low_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_amp_low (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_LG_order_Callback(hObject, eventdata, handles)
% hObject    handle to edit_LG_order (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_LG_order as text
%        str2double(get(hObject,'String')) returns contents of edit_LG_order as a double
Decon_param = getappdata(handles.figure_main,'Decon_param');
LG_order_string = get(handles.edit_LG_order,'String');
Decon_param.LG_param.order = str2double(LG_order_string);
setappdata(handles.figure_main,'Decon_param',Decon_param);


% --- Executes during object creation, after setting all properties.
function edit_LG_order_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_LG_order (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_amp_high_Callback(hObject, eventdata, handles)
% hObject    handle to edit_amp_high (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_amp_high as text
%        str2double(get(hObject,'String')) returns contents of edit_amp_high as a double
Decon_param = getappdata(handles.figure_main,'Decon_param');
amp_high_string = get(handles.edit_amp_high,'String');
Decon_param.amplitude_window(2) = str2double(amp_high_string);
setappdata(handles.figure_main,'Decon_param',Decon_param);

% --- Executes during object creation, after setting all properties.
function edit_amp_high_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_amp_high (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu_num_exp.
function popupmenu_num_exp_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_num_exp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_num_exp contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_num_exp
Decon_param = getappdata(handles.figure_main,'Decon_param');
Decon_param.Exp_param.model = get(handles.popupmenu_num_exp,'Value');
setappdata(handles.figure_main,'Decon_param',Decon_param);

% --- Executes during object creation, after setting all properties.
function popupmenu_num_exp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_num_exp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_daq_time_Callback(hObject, eventdata, handles)
% hObject    handle to edit_daq_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_daq_time as text
%        str2double(get(hObject,'String')) returns contents of edit_daq_time as a double
Decon_param = getappdata(handles.figure_main,'Decon_param');
daq_time_string = get(handles.edit_daq_time,'String');
Decon_param.DAQ_time = str2double(daq_time_string);
setappdata(handles.figure_main,'Decon_param',Decon_param);

% --- Executes during object creation, after setting all properties.
function edit_daq_time_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_daq_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu_decon_method.
function popupmenu_decon_method_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_decon_method (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_decon_method contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_decon_method
Decon_param = getappdata(handles.figure_main,'Decon_param');
Decon_param.Decon_option = get(handles.popupmenu_decon_method,'Value');
setappdata(handles.figure_main,'Decon_param',Decon_param);

% --- Executes during object creation, after setting all properties.
function popupmenu_decon_method_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_decon_method (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_ch4.
function checkbox_ch4_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_ch4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_ch4
Decon_param = getappdata(handles.figure_main,'Decon_param');
Decon_param.Channels(4) = get(handles.checkbox_ch4,'Value');
setappdata(handles.figure_main,'Decon_param',Decon_param);

% --- Executes on button press in checkbox_ch3.
function checkbox_ch3_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_ch3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_ch3
Decon_param = getappdata(handles.figure_main,'Decon_param');
Decon_param.Channels(3) = get(handles.checkbox_ch3,'Value');
setappdata(handles.figure_main,'Decon_param',Decon_param);

% --- Executes on button press in checkbox_ch2.
function checkbox_ch2_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_ch2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_ch2
Decon_param = getappdata(handles.figure_main,'Decon_param');
Decon_param.Channels(2) = get(handles.checkbox_ch2,'Value');
setappdata(handles.figure_main,'Decon_param',Decon_param);

% --- Executes on button press in checkbox_ch1.
function checkbox_ch1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_ch1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_ch1
Decon_param = getappdata(handles.figure_main,'Decon_param');
Decon_param.Channels(1) = get(handles.checkbox_ch1,'Value');
setappdata(handles.figure_main,'Decon_param',Decon_param);


% --- Executes when entered data in editable cell(s) in uitable_files.
function uitable_files_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to uitable_files (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)
Files = getappdata(handles.figure_main,'Files');
data = get(handles.uitable_files,'Data');
Files.selector = data(:,1);
setappdata(handles.figure_main,'Files',Files);
