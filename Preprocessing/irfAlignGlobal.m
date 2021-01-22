function [handles,globalShift] = irfAlignGlobal(Decon_param,name_cell,folder_cell,Background,iIRF,handles)
% tic
numOfFiles = length(name_cell);

alignInterval = 1;
ShiftAll = [];

for i=1:alignInterval:numOfFiles
%     i
    shiftTemp = zeros(1,4);
    rawDataClass = RawDataStrcutClass(name_cell(i),folder_cell(i),Background.data,iIRF.data);
    
    DC_range = getappdata(handles.figure_main,'DC_range');
    [~] = removeDC(rawDataClass,DC_range);
    
    [~, gain_list] = gainCorrection(rawDataClass);
    
    [~,~,~]= detectSaturation(rawDataClass,Decon_param.amplitude_window(1),Decon_param.amplitude_window(2));
    noise = get_noise(rawDataClass);
    if ~isempty(Background.refRange)
        [~, ~] = removeBG(rawDataClass,Background.refRange(1),Background.refRange(2));
        
    end
    peak_position = iIRF.peaks.*Decon_param.Channels;
    datastruct = truncateData(rawDataClass,Decon_param.ch_width,peak_position); %cell array to store truncated data
    
    for j=1:4
%         j
        if ~isempty(datastruct{j}')
        rawdatastruct = channeldata(datastruct{j}',iIRF.data{j},Decon_param.time_res,Decon_param.bw, '', noise,gain_list);
        Laguerre_Struct=LaguerreModel(rawdatastruct);
        Laguerre_Struct.iIRF_align();
        shiftTemp(j) = Laguerre_Struct.shift;
        end
    end
    ShiftAll = [ShiftAll;shiftTemp];
end
% shiftStd = std(ShiftAll,1);

% update main GUI field
% globalShift = mode(ShiftAll,1);
globalShift = ShiftAll;
set(handles.ch1Shift, 'String', num2str(mode(globalShift(:,1))));
set(handles.ch2Shift, 'String', num2str(mode(globalShift(:,2))));
set(handles.ch3Shift, 'String', num2str(mode(globalShift(:,3))));
set(handles.ch4Shift, 'String', num2str(mode(globalShift(:,4))));
% toc
handles.irfAlignFlag = 1; % iRF already aligned
end