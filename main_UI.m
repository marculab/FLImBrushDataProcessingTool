function varargout = main_UI(varargin)
% MAIN_UI MATLAB code for main_UI.fig
%      MAIN_UI, by itself, creates a new MAIN_UI or raises the existing
%      singleton*.
%
%      H = MAIN_UI returns the handle to a new MAIN_UI or the handle to
%      the existing singleton*.
%
%      MAIN_UI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MAIN_UI.M with the given input arguments.
%
%      MAIN_UI('Property','Value',...) creates a new MAIN_UI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before main_UI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to main_UI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help main_UI

% Last Modified by GUIDE v2.5 11-Dec-2017 19:27:46

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @main_UI_OpeningFcn, ...
    'gui_OutputFcn',  @main_UI_OutputFcn, ...
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

% ---------- Data structure attached to main_UI figure handle -----------------%
% ------------------------------------------------------------------------%
%
% Files
%   namelist
%   folderlist
%   selector
% Deconvolution parameters
%   time resolution
%   interpolation
%   bandwidth
%   channel width
%   amplitude window (low & high)
%   DAQ time ?
%   Laguerre parameters
%       Orders
%       Alpha
%   Multi-exponential parameters
%       number of exponentials
%   Phasor
%       Orders
% handles_decon
%   ~ all object handles of the Decon_UI figure
% all_data
% iIRF
%   name
%   data
%   peaks
% Background
%   reference range
%   data
% DC_range
%
%
% ------------------------------------------------------------------------%
% ---------- Data structure attached to UI figure handle -----------------%


% --- Executes just before main_UI is made visible.
function main_UI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to main_UI (see VARARGIN)

% Choose default command line output for main_UI
% main_position = get(handles.figure_main,'Position' );
% main_position(1) = 10;
% set(handles.figure_main, 'Position', main_position );

handles.output = hObject;
[main_path,~,~] = fileparts(which('main_UI.m'));
handles.main_path = main_path;
handles.rootFolder = pwd;
addpath(genpath(handles.rootFolder));
% initialize variables, change the default value here
Files = struct('folderlist',{},'namelist',{},'selector',{});
Decon_param = struct('time_res',0.08,...
    'interpolation',0.08,...
    'bw',1.5,...
    'ch_width',680,...
    'amplitude_window',[0.03 0.74],...
    'DAQ_time',9,...
    'Channels',logical([1 1 1 1]),...
    'LG_param',struct('order',12,'alpha',0.916),...
    'Exp_param',struct('model',1),... % using structure for future expandability
    'Phasor',struct('order',1,'frequency',0),... % order and frequency reserved for future development
    'Decon_option',1 ... % 1 - LG, 2 - LMEA, 3 - Exp, 4 - Phasor
    );
iIRF = struct('name','','data',[],'peaks',[]);
Background = struct('refRange',[1 2],'data',[],'dcRange',[],'path',[]);
Plot_option = struct('dataSource',1,'logScale',0,...
    'lineIndex',0,'lineMin',0,'lineMax',1,'lineStep',[1 1],...
    'pointIndex',0,'pointMin',0,'pointMax',1,'pointStep',[1 1] ...
    );
Image_option = struct('imageType',[1 2], ...
    'aspectRatio',[1 1], ...
    'pixelShift',0, ...
    'Sscan',0 ...
    );
% DC_range = [1 100];
setappdata(handles.figure_main,'Files',Files);
setappdata(handles.figure_main,'Decon_param',Decon_param);
setappdata(handles.figure_main,'irf',iIRF);
setappdata(handles.figure_main,'Background',Background);
setappdata(handles.figure_main,'Plot_option',Plot_option);
setappdata(handles.figure_main,'Image_option',Image_option);
% Launch Decon_UI module
handles_decon = Decon_UI(handles.figure_main,handles.View_decon); % start Decon_UI
setappdata(handles.figure_main,'handles_decon',handles_decon);
set(handles.View_decon,'Checked','On');
% Launch plot_UI module
handles_plot = plot_UI(handles.figure_main,handles.View_plot); % start plot UI
setappdata(handles.figure_main,'handles_plot',handles_plot);
set(handles.View_plot,'Checked','On');
% Not launching Image_UI module
setappdata(handles.figure_main,'handles_image',{});
set(handles.View_image,'Checked','Off');
handles.irfAlignFlag = 0;
%adjust location of figures
movegui(handles.figure_main,'northwest')

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes main_UI wait for user response (see UIRESUME)
% uiwait(handles.figure_main);


% --- Outputs from this function are returned to the command line.
function varargout = main_UI_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes during object deletion, before destroying properties.
function figure_main_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to figure_main (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles_decon = getappdata(handles.figure_main,'handles_decon');
if ~isempty(handles_decon)
    delete(handles_decon.figure_decon);
end
handles_plot = getappdata(handles.figure_main,'handles_plot');
if ~isempty(handles_plot)
    delete(handles_plot.figure_plot);
end
handles_image = getappdata(handles.figure_main,'handles_image');
if ~isempty(handles_image)
    delete(handles_image.figure_image);
end
close all


%**** menu Level 1 callbacks
% --------------------------------------------------------------------
function menu_file_Callback(hObject, eventdata, handles)
% hObject    handle to menu_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_irf_Callback(hObject, eventdata, handles)
% hObject    handle to menu_irf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_bg_Callback(hObject, eventdata, handles)
% hObject    handle to menu_bg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function menu_view_Callback(hObject, eventdata, handles)
% hObject    handle to menu_view (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function menu_aiming_beam_Callback(hObject, eventdata, handles)
% hObject    handle to menu_aiming_beam (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_export_Callback(hObject, eventdata, handles)
% hObject    handle to menu_export (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%**** menu Level 1 callbacks END


%**** menu Level 2 callbacks FILE
% --------------------------------------------------------------------
function file_load_Callback(hObject, eventdata, handles)
% hObject    handle to file_load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.irfAlignFlag = 0;
[namelist, folderlist] = CreatInputFileList(handles.rootFolder);
if ~isempty(namelist)
    dataRootFolder = folderlist{1};
    setappdata(handles.figure_main,'dataRootFolder',dataRootFolder);
    number_of_lines = size(namelist,1);
    selector = num2cell(logical(ones(size(namelist))));
    Files = struct('folderlist',{folderlist},'namelist',{namelist},'selector',{selector});
    setappdata(handles.figure_main,'Files',Files);
    handles_decon = getappdata(handles.figure_main,'handles_decon');
    if ~isempty(handles_decon)
        set(handles_decon.uitable_files,'Data',[Files.selector Files.namelist]);
    end
end
guidata(hObject, handles);

% --------------------------------------------------------------------
function file_load_processed_Callback(hObject, eventdata, handles)
% hObject    handle to file_load_processed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = load_processed_data(handles);

Plot_option = getappdata(handles.figure_main,'Plot_option');

all_data = getappdata(handles.figure_main,'all_data');
num_of_lines = size(all_data,1);
INT_map = getappdata(handles.figure_main,'INT_map');
num_of_point = zeros(1,4);
for kk=1:4
    num_of_point(kk) = size(INT_map{kk},2);
end
num_of_point = max(num_of_point);

if eq(num_of_lines,1)
    Plot_option.lineIndex = 0;
    Plot_option.lineMin = 0;
    Plot_option.lineMax = 1;
    Plot_option.lineStep = [1 1];
else
    Plot_option.lineIndex = 1;
    Plot_option.lineMin = 1;
    Plot_option.lineMax = num_of_lines;
    Plot_option.lineStep = [1/(num_of_lines-1) 1/(num_of_lines-1)];
end

if eq(num_of_point,1)
    Plot_option.pointIndex = 0;
    Plot_option.pointMin = 0;
    Plot_option.pointMax = 1;
    Plot_option.pointStep = [1 1];
else
    Plot_option.pointIndex = 1;
    Plot_option.pointMin = 1;
    Plot_option.pointMax = num_of_point;
    Plot_option.pointStep = [1/(num_of_point-1) 1/(num_of_point-1)];
end
setappdata(handles.figure_main,'Plot_option',Plot_option);
handles_plot = plot_UI(handles.figure_main,handles.View_plot);
setappdata(handles.figure_main,'handles_plot',handles_plot);
set(handles.View_plot,'Checked','On');
handles_image = Image_UI(handles.figure_main,handles.View_image);
setappdata(handles.figure_main,'handles_image',handles_image);
set(handles.View_image,'Checked','On');


guidata(hObject, handles);


% --------------------------------------------------------------------
function file_setpath_Callback(hObject, eventdata, handles)
% hObject    handle to file_setpath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
choice = questdlg('Which path do you want to use as root?','root path selection','GUI','Current','Select','Select');
switch choice
    case 'GUI'
        handles.rootFolder = handles.main_path;
    case 'Current'
        handles.rootFolder = pwd;
    case 'Select'
        newpath = uigetdir(handles.rootFolder,'Select new root folder');
        if ~isempty(newpath)
            handles.rootFolder = newpath;
        else
            warning('No path was selected!');
        end
    otherwise
        warning('No path was selected!');
end
% Update handles structure
guidata(hObject, handles);

%**** menu Level 2 callbacks iIRF
% --------------------------------------------------------------------
function irf_loadnew_Callback(hObject, eventdata, handles)
% hObject    handle to irf_loadnew (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.irfAlignFlag = 0;
[irf,irf_name] = loadTRIPLEXiIRF(handles.rootFolder);
iIRF = struct('name',{irf_name},'data',{irf});
setappdata(handles.figure_main,'irf',iIRF);
handles_decon = getappdata(handles.figure_main,'handles_decon');
if ~isempty(handles_decon)
    set(handles_decon.edit_iIRF_filename,'String',iIRF.name);
end

guidata(hObject, handles);

% --------------------------------------------------------------------
function irf_view_Callback(hObject, eventdata, handles)
% hObject    handle to irf_view (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
h = figure('Name','iRF window');
iIRF = getappdata(handles.figure_main,'irf');
plot(iIRF.data{1})
hold on
plot(circshift(iIRF.data{2},50));
plot(circshift(iIRF.data{3},100))
plot(circshift(iIRF.data{4},150))
hold off
title(strrep(iIRF.name,'_',' '))
legend('CH1','CH2','CH3','CH4')


%**** menu Level 2 callbacks BACKGROUND
% --------------------------------------------------------------------
function bg_loadnew_Callback(hObject, eventdata, handles)
% hObject    handle to bg_loadnew (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
dataRootFolder = getappdata(handles.figure_main,'dataRootFolder');
[BG, BG_file_path] = loadTRIPLEXBG(dataRootFolder);
BG = -BG;
handles_plot = getappdata(handles.figure_main,'handles_plot');
if ~isempty(handles_plot)
    axes(handles_plot.axes_decay)
    plot(BG);
end
% if BG_remove position is empty then no BG removal.
Background = getappdata(handles.figure_main,'Background');
if isempty(Background)
    Background = struct('refRange',[0 0],'data',BG);
else
    Background.data = BG;
    Background.path = BG_file_path;
end
setappdata(handles.figure_main,'Background',Background);

guidata(hObject, handles);

% --------------------------------------------------------------------
function bg_preprocess_Callback(hObject, eventdata, handles)
% hObject    handle to bg_preprocess (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% warning('off','all')
Decon_param= getappdata(handles.figure_main,'Decon_param');
Background = getappdata(handles.figure_main,'Background');

% get namelist and folderlist from handle
Files = getappdata(handles.figure_main,'Files');
iIRF = getappdata(handles.figure_main,'irf');
DC_range = getappdata(handles.figure_main,'DC_range');
BG_remove_flag=1;

idx = ceil(size(Files.namelist,1)/2); % chose one file from middle

temp_data_struct = RawDataStrcutClass(Files.namelist(idx),Files.folderlist(idx),Background.data,[]);
temp = get(temp_data_struct,'rawdata');

if isempty(Background.data)
    BG_skip = questdlg('No Background is loaded, do you want to skip?','BG Dialog','Load BG','Skip','Skip');
    if strcmp(BG_skip,'Load BG') %opte to load BG
        BG = -loadTRIPLEXBG(handles.rootFolder);
        Background.data = BG;
        setappdata(handles.figure_main,'Background',Background);
    else % if user chose not to load BG, display data for peak selection
        Background.data = zeros(1,size(temp,2))';
        idx = round(size(Files.namelist,1)/2);
        if idx<1
            idx=1;
        end
        if idx>size(Files.namelist,1)
            idx = size(Files.namelist,1);
        end
        
        temp_data_struct = RawDataStrcutClass(Files.namelist(idx),Files.folderlist(idx),Background.data,[]);
        temp = get(temp_data_struct,'rawdata');
        Background.data = zeros(1,size(temp,2))';
        waveToPlot = 250;
        if size(temp,1)>waveToPlot
            interval = ceil(size(temp,1)/waveToPlot);
            temp = temp(1:interval:size(temp,1),:);
        end
        
        BG_remove_fig = figure;
        set(BG_remove_fig, 'Position',  get(0,'Screensize'));
        BGPos = load(fullfile(handles.main_path,'PeaksPos.dat'));
        plot(temp');
        title('Please use green cursor for DC removal and red cursor for BG removal');
        hcur1 = mycursors(gca,'g'); % add DC subtrction cursors
        hcur1.off(1);
        hcur1.add(BGPos(4,1));
        hcur1.add(BGPos(4,2));
        
        hcur2 = mycursors(gca,'r'); % add BG removal cursors
        hcur2.off(1);
        hcur2.add(BGPos(2,1));
        hcur2.add(BGPos(2,2));
        hcur2.add(BGPos(2,3));
        hcur2.add(BGPos(2,4));
        
        button=uicontrol(BG_remove_fig,'unit','Normalized','Position',[0.425 0.09 0.15 0.05],'String','Comfirm Peak','Callback','uiresume(gcbf)');
        uiwait(BG_remove_fig);
        peak_position = [round(hcur2.val(1)) round(hcur2.val(2)) round(hcur2.val(3)) round(hcur2.val(4))];
        BGPos(2,1:4)=peak_position;
        iIRF.peaks = peak_position;
        setappdata(handles.figure_main,'irf',iIRF);
        Background.DC_range =  [round(hcur1.val(1)) round(hcur1.val(2))];
        BGPos(4,1:2) = Background.DC_range;
        setappdata(handles.figure_main,'DC_range',Background.DC_range);
        BG_remove_flag = 0;
        close(BG_remove_fig)
        save(fullfile(handles.main_path,'PeaksPos.dat'), 'BGPos', '-ascii');
    end
end

idx = round(size(Files.namelist,1)/2);
if idx<1
    idx=1;
end
if idx>size(Files.namelist,1)
    idx = size(Files.namelist,1);
end

temp_data_struct = RawDataStrcutClass(Files.namelist(idx),Files.folderlist(idx),Background.data,'');
detectSaturation(temp_data_struct,Decon_param.amplitude_window(1),Decon_param.amplitude_window(2));

idx1 = floor(temp_data_struct.rawdata_dimention(1)/2);
if idx1<1
    idx1=1;
end

temp = get(temp_data_struct,'data'); % all data
% temp_data = temp(idx1,:); % one waveform from middle of scan
waveToPlot = 250;
if size(temp,1)>waveToPlot
    interval = ceil(size(temp,1)/waveToPlot);
    temp = temp(1:interval:size(temp,1),:);
else
    interval = 1;
end

if BG_remove_flag
    fig = figure;
    set(fig, 'Position',  get(0,'Screensize'));
    % while handles.BG_remove_stop_flag
    BGPos = load(fullfile(handles.main_path,'PeaksPos.dat'));
    % plot background
    hax1 = subplot(3, 1, 1);
    
    pos = get(hax1,'position');
    pos(2) = pos(2) + 0.05;
    set(hax1,'position',pos);
    BG = Background.data;
    plot(BG);
    title('Please use GREEN cursor for DC removal and RED cursor for BG removal');
    ylim([-0.1 0.75])
    xlim([1 length(BG)]);
    hcur1 = mycursors(hax1,'g'); % add DC subtrction cursors
    hcur1.off(1);
    hcur1.add(BGPos(4,1));
    hcur1.add(BGPos(4,2));
    hcur2 = mycursors(hax1,'r');
    hcur2.off(1);
    hcur2.add(BGPos(3,1));
    hcur2.add(BGPos(3,2));
    
    % plot raw waveform
    hax2 = subplot(3, 1, 2);
    pos = get(hax2,'position');
    pos(2) = pos(2) + 0.1;
    set(hax2,'position',pos);
    % plot all data instead of only one
    plot(temp');
    xlim([1 size(temp,2)]);
    pause(0.01);
    
    button=uicontrol(fig,'unit','Normalized','Position',[0.425 0.09 0.15 0.05],'String','Remove BG','Callback','uiresume(gcbf)');
    uiwait(fig);
    Background.refRange = [ceil(min([hcur2.val(1) hcur2.val(2)])) floor(max([hcur2.val(1) hcur2.val(2)]))];
    setappdata(handles.figure_main,'Background',Background); % store BG remove cursor position
    BGPos(3,1:2) = Background.refRange;
    Background.DC_range =  [round(hcur1.val(1)) round(hcur1.val(2))];
    BGPos(4,1:2) = Background.DC_range;
    setappdata(handles.figure_main,'DC_range',Background.DC_range);
    % remove DC
    [~] = removeDC(temp_data_struct,Background.DC_range);
    plot(hax1,get(temp_data_struct,'bg'))
    ylim([-0.1 0.75])
    xlim([1 length(BG)]);
    [BG_removed, ~] = removeBG(temp_data_struct,Background.refRange(1),Background.refRange(2));
    BG_removed = BG_removed(1:interval:size(BG_removed,1),:);
    % plot BG removed waveform
    hax3 = subplot(3, 1, 3);
    pos = get(hax3,'position');
    pos(2) = pos(2) + 0.1;
    set(hax3,'position',pos);
    plot(BG_removed');
    xlim([1 size(BG_removed,2)]);
    
    hcur3 = mycursors(hax3,'r');
    hcur3.off(1);
    for ii=1:4
        if BGPos(2,ii)>size(BG_removed,2)
            BGPos(2,ii) = size(BG_removed,2)-10;
        end
    end
    hcur3.add(BGPos(2,1));
    hcur3.add(BGPos(2,2));
    hcur3.add(BGPos(2,3));
    hcur3.add(BGPos(2,4));
    
    set(button,'String','Confirm Peak');
    
    uiwait(fig)
    
    iIRF.peaks = [round(hcur3.val(1)) round(hcur3.val(2)) round(hcur3.val(3)) round(hcur3.val(4))];
    
    BGPos(2,1:4) = iIRF.peaks;
    
    setappdata(handles.figure_main,'irf',iIRF);
    close(fig)
    save(fullfile(handles.main_path,'PeaksPos.dat'), 'BGPos', '-ascii');
end
setappdata(handles.figure_main,'Background',Background);
guidata(hObject, handles);

%**** menu Level 2 callbacks VIEW
% --------------------------------------------------------------------
function View_decon_Callback(hObject, eventdata, handles)
% hObject    handle to View_decon (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles_decon = getappdata(handles.figure_main,'handles_decon');
if isempty(handles_decon)
    % Launch Decon_UI module
    handles_decon = Decon_UI(handles.figure_main,handles.View_decon);
    setappdata(handles.figure_main,'handles_decon',handles_decon);
    set(handles.View_decon,'Checked','On');
else
    delete(handles_decon.figure_decon);
end
guidata(hObject, handles);
% --------------------------------------------------------------------
function View_plot_Callback(hObject, eventdata, handles)
% hObject    handle to View_plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles_plot = getappdata(handles.figure_main,'handles_plot');
if isempty(handles_plot)
    % Launch plot_UI module
    handles_plot = plot_UI(handles.figure_main,handles.View_plot);
    setappdata(handles.figure_main,'handles_plot',handles_plot);
    set(handles.View_plot,'Checked','On');
else
    delete(handles_plot.figure_plot);
end
guidata(hObject, handles);
% --------------------------------------------------------------------
function View_image_Callback(hObject, eventdata, handles)
% hObject    handle to View_image (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles_image = getappdata(handles.figure_main,'handles_image');
if isempty(handles_image)
    % Launch Image_UI module
    handles_image = Image_UI(handles.figure_main,handles.View_image);
    setappdata(handles.figure_main,'handles_image',handles_image);
    set(handles.View_image,'Checked','On');
else
    delete(handles_image.figure_image);
end
guidata(hObject, handles);

%**** menu Level 2 callbacks AIMING BEAM
% --------------------------------------------------------------------
function ab_load_Callback(hObject, eventdata, handles)
% hObject    handle to ab_load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function ab_segmentation_Callback(hObject, eventdata, handles)
% hObject    handle to ab_segmentation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


%**** menu Level 2 callbacks EXPORT
% --------------------------------------------------------------------
function export_csv_Callback(hObject, eventdata, handles)
% hObject    handle to export_csv (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


%**** menu Level 2 callbacks END

%**** button callbacks
% --- Executes on button press in btn_proc_all.
function btn_proc_all_Callback(hObject, eventdata, handles)
% hObject    handle to btn_proc_all (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
warning('off')
set(handles.ch1Shift, 'enable', 'off');
set(handles.ch2Shift, 'enable', 'off');
set(handles.ch3Shift, 'enable', 'off');
set(handles.ch4Shift, 'enable', 'off');
% get namelist and folderlist from handle
Files = getappdata(handles.figure_main,'Files');
% ger IRF
iIRF = getappdata(handles.figure_main,'irf');
Decon_param = getappdata(handles.figure_main,'Decon_param');
Background = getappdata(handles.figure_main,'Background');
Plot_option = getappdata(handles.figure_main,'Plot_option');

if ~isempty(iIRF.peaks)
    if ~any(iIRF.peaks)
        msgbox('No data will be deconvolved, please choose the righr channel or reselect peak position.');
    end
else
    msgbox('No peak position is selected, please do Background subtraction and peak selection.')
    error('No peak_position')
end

[name_cell, folder_cell] = creatJobQueue(Files.namelist, Files.folderlist,1);
handles.jPb = javax.swing.JProgressBar;
set(handles.jPb,'StringPainted',1,'Value',12.5,'Indeterminate',0);
[handles.hPb, handles.hContainer] = javacomponent(handles.jPb,[550 45 200 20],handles.figure_main);

switch Decon_param.Decon_option
    case 1
        handles = run_Deconvolution(Decon_param,name_cell,folder_cell,Background,iIRF,handles,0);
    case 3
        handles = run_Deconvolution(Decon_param,name_cell,folder_cell,Background,iIRF,handles,0);
    case 2
        handles = run_Decon(Decon_param,Files,Background,iIRF,handles);
    case 5
        handles = run_Deconvolution(Decon_param,name_cell,folder_cell,Background,iIRF,handles,0);
    otherwise
        error('Unknown deconvolution method!');
end

if handles.decon_success_flag
    set(handles.edit_decon_process,'String','Finished');
    drawnow
    
    all_data = getappdata(handles.figure_main,'all_data');
    num_of_lines = size(all_data,1);
    INT_map = getappdata(handles.figure_main,'INT_map');
    num_of_point = zeros(1,4);
    for kk=1:4
        num_of_point(kk) = size(INT_map{kk},2);
    end
    num_of_point = max(num_of_point);
    
    if eq(num_of_lines,1)
        Plot_option.lineIndex = 0;
        Plot_option.lineMin = 0;
        Plot_option.lineMax = 1;
        Plot_option.lineStep = [1 1];
    else
        Plot_option.lineIndex = 1;
        Plot_option.lineMin = 1;
        Plot_option.lineMax = num_of_lines;
        Plot_option.lineStep = [1/(num_of_lines-1) 1/(num_of_lines-1)];
    end
    
    if eq(num_of_point,1)
        Plot_option.pointIndex = 0;
        Plot_option.pointMin = 0;
        Plot_option.pointMax = 1;
        Plot_option.pointStep = [1 1];
    else
        Plot_option.pointIndex = 1;
        Plot_option.pointMin = 1;
        Plot_option.pointMax = num_of_point;
        Plot_option.pointStep = [1/(num_of_point-1) 1/(num_of_point-1)];
    end
    setappdata(handles.figure_main,'Plot_option',Plot_option);
    handles_plot = plot_UI(handles.figure_main,handles.View_plot);
    setappdata(handles.figure_main,'handles_plot',handles_plot);
    set(handles.View_plot,'Checked','On');
    handles_image = Image_UI(handles.figure_main,handles.View_image);
    setappdata(handles.figure_main,'handles_image',handles_image);
    set(handles.View_image,'Checked','On');
end

set(handles.ch1Shift, 'enable', 'on');
set(handles.ch2Shift, 'enable', 'on');
set(handles.ch3Shift, 'enable', 'on');
set(handles.ch4Shift, 'enable', 'on');
guidata(hObject, handles);


% --- Executes on button press in btn_proc_select.
function btn_proc_select_Callback(hObject, eventdata, handles)
% hObject    handle to btn_proc_select (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Files = getappdata(handles.figure_main,'Files');
if ~isempty(Files)
    warning('off')
    iIRF = getappdata(handles.figure_main,'irf');
    Decon_param = getappdata(handles.figure_main,'Decon_param');
    Background = getappdata(handles.figure_main,'Background');
    Plot_option = getappdata(handles.figure_main,'Plot_option');
    
    if ~isempty(iIRF.peaks)
        if ~any(iIRF.peaks)
            msgbox('No data will be deconvolved, please choose the righr channel or reselect peak position.');
        end
    else
        msgbox('No peak position is selected, please do Background subtraction and peak selection.')
        error('No peak_position')
    end
    
    [name_cell, folder_cell] = creatJobQueue(Files.namelist(cell2mat(Files.selector)), Files.folderlist(cell2mat(Files.selector)),1);
    handles.jPb = javax.swing.JProgressBar;
    set(handles.jPb,'StringPainted',1,'Value',12.5,'Indeterminate',0);
    
    [handles.hPb, handles.hContainer] = javacomponent(handles.jPb,[550 5 200 20],handles.figure_main);
    
    handles = run_Deconvolution(Decon_param,name_cell,folder_cell,Background,iIRF,handles,0);
    
    
    if handles.decon_success_flag
        set(handles.edit_decon_process,'String','Finished');
        drawnow
        
        all_data = getappdata(handles.figure_main,'all_data');
        num_of_lines = size(all_data,1);
        INT_map = getappdata(handles.figure_main,'INT_map');
        num_of_point = zeros(1,4);
        for kk=1:4
            num_of_point(kk) = size(INT_map{kk},2);
        end
        num_of_point = max(num_of_point);
        if eq(num_of_lines,1)
            Plot_option.lineIndex = 0;
            Plot_option.lineMin = 0;
            Plot_option.lineMax = 1;
            Plot_option.lineStep = [1 1];
        else
            Plot_option.lineIndex = 1;
            Plot_option.lineMin = 1;
            Plot_option.lineMax = num_of_lines;
            Plot_option.lineStep = [1/(num_of_lines-1) 1/(num_of_lines-1)];
        end
        
        if eq(num_of_point,1)
            Plot_option.pointIndex = 0;
            Plot_option.pointMin = 0;
            Plot_option.pointMax = 1;
            Plot_option.pointStep = [1 1];
        else
            Plot_option.pointIndex = 1;
            Plot_option.pointMin = 1;
            Plot_option.pointMax = num_of_point;
            Plot_option.pointStep = [1/(num_of_point-1) 1/(num_of_point-1)];
        end
        setappdata(handles.figure_main,'Plot_option',Plot_option);
        handles_plot = plot_UI(handles.figure_main,handles.View_plot);
        setappdata(handles.figure_main,'handles_plot',handles_plot);
        set(handles.View_plot,'Checked','On');
        handles_image = Image_UI(handles.figure_main,handles.View_image);
        setappdata(handles.figure_main,'handles_image',handles_image);
        set(handles.View_image,'Checked','On');
        
    end % END OF ~isempty(Files)
end
guidata(hObject, handles);

% --- Executes on button press in btn_stop.
function btn_stop_Callback(hObject, eventdata, handles)
% hObject    handle to btn_stop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% the Stop button is reserved for parallel processing version

%**** button callbacks END


% --------------------------------------------------------------------
function export_mat_Callback(hObject, eventdata, handles)
% hObject    handle to export_mat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
warning('off')
% get namelist and folderlist from handle
Files = getappdata(handles.figure_main,'Files');
% ger IRF
iIRF = getappdata(handles.figure_main,'irf');
Decon_param = getappdata(handles.figure_main,'Decon_param');
Background = getappdata(handles.figure_main,'Background');
Plot_option = getappdata(handles.figure_main,'Plot_option');

if ~isempty(iIRF.peaks)
    if ~any(iIRF.peaks)
        msgbox('No data will be deconvolved, please choose the righr channel or reselect peak position.');
    end
else
    msgbox('No peak position is selected, please do Background subtraction and peak selection.')
    error('No peak_position')
end

[name_cell, folder_cell] = creatJobQueue(Files.namelist, Files.folderlist,1);
handles.jPb = javax.swing.JProgressBar;
set(handles.jPb,'StringPainted',1,'Value',12.5,'Indeterminate',0);
[handles.hPb, handles.hContainer] = javacomponent(handles.jPb,[550 5 200 20],handles.figure_main);

handles = exportPreprossedDataNew(Decon_param,name_cell,folder_cell,Background,iIRF,handles,0);

guidata(hObject, handles);

% --------------------------------------------------------------------
function export_tiff_Callback(hObject, eventdata, handles)
% hObject    handle to export_tiff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

warning('off')
% get namelist and folderlist from handle
Files = getappdata(handles.figure_main,'Files');
% ger IRF
iIRF = getappdata(handles.figure_main,'irf');
Decon_param = getappdata(handles.figure_main,'Decon_param');
Background = getappdata(handles.figure_main,'Background');
Plot_option = getappdata(handles.figure_main,'Plot_option');

if ~isempty(iIRF.peaks)
    if ~any(iIRF.peaks)
        msgbox('No data will be deconvolved, please choose the righr channel or reselect peak position.');
    end
else
    msgbox('No peak position is selected, please do Background subtraction and peak selection.')
    error('No peak_position')
end

[name_cell, folder_cell] = creatJobQueue(Files.namelist, Files.folderlist,1);
handles.jPb = javax.swing.JProgressBar;
set(handles.jPb,'StringPainted',1,'Value',12.5,'Indeterminate',0);
[handles.hPb, handles.hContainer] = javacomponent(handles.jPb,[550 5 200 20],handles.figure_main);

handles = exportTiffStackNew(Decon_param,name_cell,folder_cell,Background,iIRF,handles,0);
guidata(hObject, handles);

% --------------------------------------------------------------------
function export_ffh_Callback(hObject, eventdata, handles)
% hObject    handle to export_ffh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
warning('off')
% get namelist and folderlist from handle
Files = getappdata(handles.figure_main,'Files');
% ger IRF
iIRF = getappdata(handles.figure_main,'irf');
Decon_param = getappdata(handles.figure_main,'Decon_param');
Background = getappdata(handles.figure_main,'Background');
Plot_option = getappdata(handles.figure_main,'Plot_option');

if ~isempty(iIRF.peaks)
    if ~any(iIRF.peaks)
        msgbox('No data will be deconvolved, please choose the righr channel or reselect peak position.');
    end
else
    msgbox('No peak position is selected, please do Background subtraction and peak selection.')
    error('No peak_position')
end

[name_cell, folder_cell] = creatJobQueue(Files.namelist, Files.folderlist,1);
handles.jPb = javax.swing.JProgressBar;
set(handles.jPb,'StringPainted',1,'Value',12.5,'Indeterminate',0);
[handles.hPb, handles.hContainer] = javacomponent(handles.jPb,[550 300 200 20],handles.figure_main);

handles = exportFLIMFitFFH(Decon_param,name_cell,folder_cell,Background,iIRF,handles,0);


% --- Executes on button press in btn_global_align.
function btn_global_align_Callback(hObject, eventdata, handles)
% hObject    handle to btn_global_align (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% handles = run_Deconvolution(Decon_param,name_cell,folder_cell,Background,iIRF,handles,0);
Files = getappdata(handles.figure_main,'Files');
% ger IRF
iIRF = getappdata(handles.figure_main,'irf');
Decon_param = getappdata(handles.figure_main,'Decon_param');
Background = getappdata(handles.figure_main,'Background');
Plot_option = getappdata(handles.figure_main,'Plot_option');

if ~isempty(iIRF.peaks)
    if ~any(iIRF.peaks)
        msgbox('No data will be deconvolved, please choose the righr channel or reselect peak position.');
    end
else
    msgbox('No peak position is selected, please do Background subtraction and peak selection.')
    error('No peak_position')
end

[name_cell, folder_cell] = creatJobQueue(Files.namelist, Files.folderlist,1);

temp_text = 'Global iRF alignment in process.';
set(handles.edit_decon_process,'String',temp_text);

drawnow

[handles, global_shift] = irfAlignGlobal(Decon_param,name_cell,folder_cell,Background,iIRF,handles);
%save iRF align result-----------------------------------------------------
cd(folder_cell{1});
mkdir('GlobaliRFAlign')
cd('GlobaliRFAlign')
fileID = fopen('iRF_Alignment.txt','w');
fprintf(fileID,'%6s %6s %6s %6s\r\n','CH1','CH2','CH3','CH4');
fprintf(fileID,'%6.0f %6.0f %6.0f %6.0f\r\n',global_shift');
fclose(fileID);
%---------------------------------------------------------------------------
set(handles.ch1Shift, 'enable', 'on');
set(handles.ch2Shift, 'enable', 'on');
set(handles.ch3Shift, 'enable', 'on');
set(handles.ch4Shift, 'enable', 'on');
msgbox('Global Alignment Finished!')

set(handles.ch1Shift, 'enable', 'on');
set(handles.ch2Shift, 'enable', 'on');
set(handles.ch3Shift, 'enable', 'on');
set(handles.ch4Shift, 'enable', 'on');

% handles.irfAlignFlag = 1; % iRF already aligned
guidata(hObject, handles);



function ch1Shift_Callback(hObject, eventdata, handles)
% hObject    handle to ch1Shift (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ch1Shift as text
%        str2double(get(hObject,'String')) returns contents of ch1Shift as a double


% --- Executes during object creation, after setting all properties.
function ch1Shift_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ch1Shift (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ch2Shift_Callback(hObject, eventdata, handles)
% hObject    handle to ch2Shift (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ch2Shift as text
%        str2double(get(hObject,'String')) returns contents of ch2Shift as a double


% --- Executes during object creation, after setting all properties.
function ch2Shift_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ch2Shift (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ch3Shift_Callback(hObject, eventdata, handles)
% hObject    handle to ch3Shift (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ch3Shift as text
%        str2double(get(hObject,'String')) returns contents of ch3Shift as a double


% --- Executes during object creation, after setting all properties.
function ch3Shift_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ch3Shift (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ch4Shift_Callback(hObject, eventdata, handles)
% hObject    handle to ch4Shift (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ch4Shift as text
%        str2double(get(hObject,'String')) returns contents of ch4Shift as a double


% --- Executes during object creation, after setting all properties.
function ch4Shift_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ch4Shift (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
