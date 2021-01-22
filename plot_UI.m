function varargout = plot_UI(varargin)
% PLOT_UI MATLAB code for plot_UI.fig
%      PLOT_UI, by itself, creates a new PLOT_UI or raises the existing
%      singleton*.
%
%      H = PLOT_UI returns the handle to a new PLOT_UI or the handle to
%      the existing singleton*.
%
%      PLOT_UI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PLOT_UI.M with the given input arguments.
%
%      PLOT_UI('Property','Value',...) creates a new PLOT_UI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before plot_UI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to plot_UI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help plot_UI

% Last Modified by GUIDE v2.5 22-May-2018 11:27:01

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @plot_UI_OpeningFcn, ...
                   'gui_OutputFcn',  @plot_UI_OutputFcn, ...
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


% --- Executes just before plot_UI is made visible.
function plot_UI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to plot_UI (see VARARGIN)

movegui(handles.figure_plot,'center');
position = get(handles.figure_plot,'Position');
position(2) = position(2)-30;
position(1) = position(1)-100;
set(handles.figure_plot, 'Position', position);

handles.figure_main = varargin{1};
handles.View_plot = varargin{2};

Plot_option = getappdata(handles.figure_main,'Plot_option');
set(handles.popupmenu_data_source,'Value',Plot_option.dataSource);
set(handles.checkbox_log_scale,'Value',Plot_option.logScale);
set(handles.slider_line,'Value',Plot_option.lineIndex,...
    'Min',Plot_option.lineMin,...
    'Max',Plot_option.lineMax,...
    'SliderStep',Plot_option.lineStep...
    );
set(handles.slider_point,'Value',Plot_option.pointIndex,...
    'Min',Plot_option.pointMin,...
    'Max',Plot_option.pointMax,...
    'SliderStep',Plot_option.pointStep...
    );

if ~isempty(getappdata(handles.figure_main,'all_data'))
    handles = Plot_DeCon_Result(handles);
end

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes plot_UI wait for user response (see UIRESUME)
% uiwait(handles.figure_plot);


% --- Outputs from this function are returned to the command line.
function varargout = plot_UI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles;


% --- Executes during object deletion, before destroying properties.
function figure_plot_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to figure_plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.View_plot,'Checked','Off');
setappdata(handles.figure_main,'handles_plot',{});


% --- Executes on selection change in popupmenu_data_source.
function popupmenu_data_source_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_data_source (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_data_source contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_data_source
Plot_option = getappdata(handles.figure_main,'Plot_option');
Plot_option.dataSource = get(handles.popupmenu_data_source,'Value');
setappdata(handles.figure_main,'Plot_option',Plot_option);

if ~isempty(getappdata(handles.figure_main,'all_data'))
    handles = Plot_DeCon_Result(handles);
end

guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function popupmenu_data_source_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_data_source (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_log_scale.
function checkbox_log_scale_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_log_scale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_log_scale
Plot_option = getappdata(handles.figure_main,'Plot_option');
Plot_option.logScale = get(handles.checkbox_log_scale,'Value');
setappdata(handles.figure_main,'Plot_option',Plot_option);

if ~isempty(getappdata(handles.figure_main,'all_data'))
    handles = Plot_DeCon_Result(handles);
end

guidata(hObject, handles);


% --- Executes on slider movement.
function slider_point_Callback(hObject, eventdata, handles)
% hObject    handle to slider_point (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
Plot_option = getappdata(handles.figure_main,'Plot_option');
Plot_option.pointIndex = floor(get(handles.slider_point,'Value'));
setappdata(handles.figure_main,'Plot_option',Plot_option);

if ~isempty(getappdata(handles.figure_main,'all_data'))
    handles = Plot_DeCon_Result(handles);
end

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function slider_point_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_point (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider_line_Callback(hObject, eventdata, handles)
% hObject    handle to slider_line (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
Plot_option = getappdata(handles.figure_main,'Plot_option');
Plot_option.lineIndex = floor(get(handles.slider_line,'Value'));
setappdata(handles.figure_main,'Plot_option',Plot_option);

if ~isempty(getappdata(handles.figure_main,'all_data'))
    handles = Plot_DeCon_Result(handles);
end

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function slider_line_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_line (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function x_max_Callback(hObject, eventdata, handles)
% hObject    handle to x_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of x_max as text
%        str2double(get(hObject,'String')) returns contents of x_max as a double
Plot_option = getappdata(handles.figure_main,'Plot_option');
Plot_option.lineIndex = floor(get(handles.slider_line,'Value'));
setappdata(handles.figure_main,'Plot_option',Plot_option);

if ~isempty(getappdata(handles.figure_main,'all_data'))
    handles = Plot_DeCon_Result(handles);
end

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function x_max_CreateFcn(hObject, eventdata, handles)
% hObject    handle to x_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function x_min_Callback(hObject, eventdata, handles)
% hObject    handle to x_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of x_min as text
%        str2double(get(hObject,'String')) returns contents of x_min as a double
Plot_option = getappdata(handles.figure_main,'Plot_option');
Plot_option.lineIndex = floor(get(handles.slider_line,'Value'));
setappdata(handles.figure_main,'Plot_option',Plot_option);

if ~isempty(getappdata(handles.figure_main,'all_data'))
    handles = Plot_DeCon_Result(handles);
end

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function x_min_CreateFcn(hObject, eventdata, handles)
% hObject    handle to x_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when figure_plot is resized.
function figure_plot_SizeChangedFcn(hObject, eventdata, handles)
% hObject    handle to figure_plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
