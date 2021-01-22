function varargout = Image_UI(varargin)
% IMAGE_UI MATLAB code for Image_UI.fig
%      IMAGE_UI, by itself, creates a new IMAGE_UI or raises the existing
%      singleton*.
%
%      H = IMAGE_UI returns the handle to a new IMAGE_UI or the handle to
%      the existing singleton*.
%
%      IMAGE_UI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in IMAGE_UI.M with the given input arguments.
%
%      IMAGE_UI('Property','Value',...) creates a new IMAGE_UI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Image_UI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Image_UI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Image_UI

% Last Modified by GUIDE v2.5 21-May-2018 12:21:57

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Image_UI_OpeningFcn, ...
                   'gui_OutputFcn',  @Image_UI_OutputFcn, ...
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


% --- Executes just before Image_UI is made visible.
function Image_UI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Image_UI (see VARARGIN)

handles.figure_main = varargin{1};
handles.View_image = varargin{2};
movegui(handles.figure_image,'southeast');
Image_option = getappdata(handles.figure_main,'Image_option');
set(handles.popupmenu_top,'Value',Image_option.imageType(1));
set(handles.popupmenu_bottom,'Value',Image_option.imageType(2));
set(handles.edit_xaspect,'String',sprintf('%g',Image_option.aspectRatio(1)));
set(handles.edit_yaspect,'String',sprintf('%g',Image_option.aspectRatio(2)));
set(handles.edit_pixel_shift,'String',sprintf('%g',Image_option.pixelShift));
set(handles.checkbox_sscan,'Value',Image_option.Sscan);

if ~isempty(getappdata(handles.figure_main,'all_data'))
    handles = plot_maps(handles,'top');
    handles = plot_maps(handles,'bottom');
end

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Image_UI wait for user response (see UIRESUME)
% uiwait(handles.figure_image);


% --- Outputs from this function are returned to the command line.
function varargout = Image_UI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles;


% --- Executes during object deletion, before destroying properties.
function figure_image_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to figure_image (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.View_image,'Checked','Off');
setappdata(handles.figure_main,'handles_image',{});

% --- Executes on selection change in popupmenu_top.
function popupmenu_top_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_top (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_top contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_top
Image_option = getappdata(handles.figure_main,'Image_option');
Image_option.imageType(1) = get(handles.popupmenu_top,'Value');
setappdata(handles.figure_main,'Image_option',Image_option);

if ~isempty(getappdata(handles.figure_main,'all_data'))
    handles = plot_maps(handles,'top');
    handles = plot_maps(handles,'bottom');
end

guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function popupmenu_top_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_top (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu_bottom.
function popupmenu_bottom_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_bottom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_bottom contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_bottom
Image_option = getappdata(handles.figure_main,'Image_option');
Image_option.imageType(2) = get(handles.popupmenu_bottom,'Value');
setappdata(handles.figure_main,'Image_option',Image_option);

if ~isempty(getappdata(handles.figure_main,'all_data'))
    handles = plot_maps(handles,'top');
    handles = plot_maps(handles,'bottom');
end

guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function popupmenu_bottom_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_bottom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_yaspect_Callback(hObject, eventdata, handles)
% hObject    handle to edit_yaspect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_yaspect as text
%        str2double(get(hObject,'String')) returns contents of edit_yaspect as a double
Image_option = getappdata(handles.figure_main,'Image_option');
Image_option.aspectRatio(2) = str2double(get(handles.edit_yaspect,'String'));
setappdata(handles.figure_main,'Image_option',Image_option);

if ~isempty(getappdata(handles.figure_main,'all_data'))
    handles = plot_maps(handles,'top');
    handles = plot_maps(handles,'bottom');
end

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit_yaspect_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_yaspect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_xaspect_Callback(hObject, eventdata, handles)
% hObject    handle to edit_xaspect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_xaspect as text
%        str2double(get(hObject,'String')) returns contents of edit_xaspect as a double
Image_option = getappdata(handles.figure_main,'Image_option');
Image_option.aspectRatio(1) = str2double(get(handles.edit_xaspect,'String'));
setappdata(handles.figure_main,'Image_option',Image_option);

if ~isempty(getappdata(handles.figure_main,'all_data'))
    handles = plot_maps(handles,'top');
    handles = plot_maps(handles,'bottom');
end

guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function edit_xaspect_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_xaspect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_pixel_shift_Callback(hObject, eventdata, handles)
% hObject    handle to edit_pixel_shift (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_pixel_shift as text
%        str2double(get(hObject,'String')) returns contents of edit_pixel_shift as a double
Image_option = getappdata(handles.figure_main,'Image_option');
Image_option.pixelShift = str2double(get(handles.edit_pixel_shift,'String'));
setappdata(handles.figure_main,'Image_option',Image_option);

if ~isempty(getappdata(handles.figure_main,'all_data'))
    handles = plot_maps(handles,'top');
    handles = plot_maps(handles,'bottom');
end

guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function edit_pixel_shift_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_pixel_shift (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_sscan.
function checkbox_sscan_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_sscan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_sscan
Image_option = getappdata(handles.figure_main,'Image_option');
Image_option.Sscan = get(handles.checkbox_sscan,'Value');
setappdata(handles.figure_main,'Image_option',Image_option);

if ~isempty(getappdata(handles.figure_main,'all_data'))
    handles = plot_maps(handles,'top');
    handles = plot_maps(handles,'bottom');
end

guidata(hObject, handles);



function SNR_Threshold_Callback(hObject, eventdata, handles)
% hObject    handle to SNR_Threshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of SNR_Threshold as text
%        str2double(get(hObject,'String')) returns contents of SNR_Threshold as a double


% --- Executes during object creation, after setting all properties.
function SNR_Threshold_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SNR_Threshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function LT_low_Callback(hObject, eventdata, handles)
% hObject    handle to LT_low (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of LT_low as text
%        str2double(get(hObject,'String')) returns contents of LT_low as a double


% --- Executes during object creation, after setting all properties.
function LT_low_CreateFcn(hObject, eventdata, handles)
% hObject    handle to LT_low (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function LT_high_Callback(hObject, eventdata, handles)
% hObject    handle to LT_high (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of LT_high as text
%        str2double(get(hObject,'String')) returns contents of LT_high as a double


% --- Executes during object creation, after setting all properties.
function LT_high_CreateFcn(hObject, eventdata, handles)
% hObject    handle to LT_high (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Interp_check_box.
function Interp_check_box_Callback(hObject, eventdata, handles)
% hObject    handle to Interp_check_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Interp_check_box
