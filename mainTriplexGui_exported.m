classdef mainTriplexGui_exported < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                        matlab.ui.Figure
        MainApp                         matlab.ui.container.Menu
        LoadRawDataMenu                 matlab.ui.container.Menu
        LoadProcessedDataMenu           matlab.ui.container.Menu
        SetRootPathMenu                 matlab.ui.container.Menu
        ExportMenu                      matlab.ui.container.Menu
        DataFileSelectedDisp            matlab.ui.control.EditField
        DataFileSelectedLabel           matlab.ui.control.Label
        TabGroup                        matlab.ui.container.TabGroup
        fileSelectionTab                matlab.ui.container.Tab
        GridLayout                      matlab.ui.container.GridLayout
        dataFilesLabel                  matlab.ui.control.Label
        dataFielesListBox               matlab.ui.control.ListBox
        apd1LoadedLamp                  matlab.ui.control.Lamp
        CalibrationLabel                matlab.ui.control.Label
        caliFileEditField               matlab.ui.control.EditField
        APD3Label                       matlab.ui.control.Label
        apd3FileEditField               matlab.ui.control.EditField
        APD2Label                       matlab.ui.control.Label
        apd2FileEditField               matlab.ui.control.EditField
        BGfileEditFieldLabel            matlab.ui.control.Label
        bgFileEditField                 matlab.ui.control.EditField
        APD1Label                       matlab.ui.control.Label
        apd1FileEditField               matlab.ui.control.EditField
        bgLoadedLamp                    matlab.ui.control.Lamp
        apd1LoadButton                  matlab.ui.control.Button
        bgLoadButton                    matlab.ui.control.Button
        bgNotNeededCheckBox             matlab.ui.control.CheckBox
        loadDataButton                  matlab.ui.control.Button
        clearButton                     matlab.ui.control.Button
        apd2LoadedLamp                  matlab.ui.control.Lamp
        apd2LoadButton                  matlab.ui.control.Button
        apd3LoadedLamp                  matlab.ui.control.Lamp
        apd3LoadButton                  matlab.ui.control.Button
        apd1ViewButton                  matlab.ui.control.Button
        apd2ViewButton                  matlab.ui.control.Button
        apd3ViewButton                  matlab.ui.control.Button
        caliLoadedLamp                  matlab.ui.control.Lamp
        caliLoadButton                  matlab.ui.control.Button
        caliViewButton                  matlab.ui.control.Button
        bg3UIAxes                       matlab.ui.control.UIAxes
        bg2UIAxes                       matlab.ui.control.UIAxes
        bg1UIAxes                       matlab.ui.control.UIAxes
        DeConSettingTab                 matlab.ui.container.Tab
        GridLayout2                     matlab.ui.container.GridLayout
        RunEXPCheckBox                  matlab.ui.control.CheckBox
        LogScaleCheckBox                matlab.ui.control.CheckBox
        TruncationPanel                 matlab.ui.container.Panel
        Ch3AlphaValueDropDown           matlab.ui.control.DropDown
        Ch3AlphaValueDropDownLabel      matlab.ui.control.Label
        Ch2AlphaValueDropDown           matlab.ui.control.DropDown
        Ch2AlphaValueDropDownLabel      matlab.ui.control.Label
        Ch1AlphaValueDropDown           matlab.ui.control.DropDown
        Ch1AlphaValueDropDownLabel      matlab.ui.control.Label
        LaguerreSettingPanel            matlab.ui.container.Panel
        exponentialsDropDown            matlab.ui.control.DropDown
        ofexponentialsDropDownLabel     matlab.ui.control.Label
        MultiExponentialSettingPanel    matlab.ui.container.Panel
        IgnoreDataRangePanel            matlab.ui.container.Panel
        BackgroudRemovalPanel           matlab.ui.container.Panel
        ChannelSelectionPanel           matlab.ui.container.Panel
        channel3CheckBox                matlab.ui.control.CheckBox
        channel2CheckBox                matlab.ui.control.CheckBox
        channel1CheckBox                matlab.ui.control.CheckBox
        DigitizerSettingPanel           matlab.ui.container.Panel
        DataAveragingPanel              matlab.ui.container.Panel
        SignalConditioningPanel         matlab.ui.container.Panel
        HighThresholdEditFieldLabel     matlab.ui.control.Label
        HighThresholdEditField          matlab.ui.control.NumericEditField
        LowThresholdEditFieldLabel      matlab.ui.control.Label
        LowThresholdEditField           matlab.ui.control.NumericEditField
        DataAverageDropDown             matlab.ui.control.DropDown
        DataAverageDropDownLabel        matlab.ui.control.Label
        TimeResolutionnsEditFieldLabel  matlab.ui.control.Label
        TimeResolutionnsEditField       matlab.ui.control.NumericEditField
        BgThreholdEditFieldLabel        matlab.ui.control.Label
        BgThreholdEditField             matlab.ui.control.NumericEditField
        BgHighEditFieldLabel            matlab.ui.control.Label
        BgHighEditField                 matlab.ui.control.NumericEditField
        BgLowEditFieldLabel             matlab.ui.control.Label
        BgLowEditField                  matlab.ui.control.NumericEditField
        DcHighEditFieldLabel            matlab.ui.control.Label
        dataIgnoreHighEditField         matlab.ui.control.NumericEditField
        DcLowEditFieldLabel             matlab.ui.control.Label
        dataIgnoreLowEditField          matlab.ui.control.NumericEditField
        LaguerreOrderEditFieldLabel     matlab.ui.control.Label
        LaguerreOrderEditField          matlab.ui.control.NumericEditField
        ChannelWidthnsLabel             matlab.ui.control.Label
        channelWidthEditField           matlab.ui.control.NumericEditField
        LoadDataButton                  matlab.ui.control.Button
        RemoveBGButton                  matlab.ui.control.Button
        UpsampleButton                  matlab.ui.control.Button
        laserRepRateEditField           matlab.ui.control.NumericEditField
        nmlaserreprateEditFieldLabel    matlab.ui.control.Label
        DeconButton                     matlab.ui.control.Button
        AverageDataButton               matlab.ui.control.Button
        TruncateButton                  matlab.ui.control.Button
        RemoveSaturationButton          matlab.ui.control.Button
        UIAxes_ch3BgRemoveFig           matlab.ui.control.UIAxes
        UIAxes_ch2BgRemoveFig           matlab.ui.control.UIAxes
        UIAxes_ch3DataFig               matlab.ui.control.UIAxes
        UIAxes_ch1BgRemoveFig           matlab.ui.control.UIAxes
        UIAxes_ch2DataFig               matlab.ui.control.UIAxes
        UIAxes_ch1DataFig               matlab.ui.control.UIAxes
        FittingResultTab                matlab.ui.container.Tab
        GridLayout3                     matlab.ui.container.GridLayout
        PointEditField                  matlab.ui.control.NumericEditField
        PointEditFieldLabel             matlab.ui.control.Label
        LineSlider                      matlab.ui.control.Slider
        PlotLabel                       matlab.ui.control.Label
        SaveDataButton                  matlab.ui.control.Button
        xMaxEditFieldLabel              matlab.ui.control.Label
        xMaxEditField                   matlab.ui.control.NumericEditField
        xMinEditFieldLabel              matlab.ui.control.Label
        xMinEditField                   matlab.ui.control.NumericEditField
        PlotDropDown                    matlab.ui.control.DropDown
        channelLabel                    matlab.ui.control.Label
        channelDropDown                 matlab.ui.control.DropDown
        UIAxesDeconResult               matlab.ui.control.UIAxes
        UIAxesResidue                   matlab.ui.control.UIAxes
        UIAxesAutoCo                    matlab.ui.control.UIAxes
        UIAxesFitting                   matlab.ui.control.UIAxes
        imagesReconstrctionTab          matlab.ui.container.Tab
        EditField                       matlab.ui.control.EditField
        EditFieldLabel                  matlab.ui.control.Label
        UIAxesImg                       matlab.ui.control.UIAxes
        ContextMenu                     matlab.ui.container.ContextMenu
    end


    properties (Access = private)
        oldpath = addpath(genpath(pwd));

        apd1DialogApp % apd1 dialog app
        apd2DialogApp % apd1 dialog app
        apd3DialogApp % apd1 dialog app

        poolobj % parallel pool object
        digitizerNoise = 0.0078; % ENOB = 8;
        EOP_H1G % Extended output phasor real
        EOP_H1S % Extended output phasor imag
        SP_G % spectral phasor real
        SP_S % spectral phasor imag
    end

    properties (Access = public)
        upSampleFactor = 5;
        irfAlignFlag = 0; % irf alignment flag
        apd1LoadedFlag =0; % apd1 load flag
        apd2LoadedFlag =0; % apd2 load flag
        apd3LoadedFlag =0; % apd3 load flag
        bgLoadedFlag = 0; % bg load flag
        caliLoadedFlag = 0; % calibratin loaded flag
        expDeconFlag = 0; % Multi-exp deconvolusion flag
        processedDataloadedFlag = 0; % Description
        rootFolder = fileparts(mfilename('fullpath')) % current working directory
        dataInfoObj = dataInfo('','','','','','','','','','','','');% Description
        dataNameList = [];
        dataFolder = pwd;
        detectorDataFolder = pwd;
        apd1Obj = [];
        apd2Obj = [];
        apd3Obj = [];
        irfName = [];
        bgObj = [];
        caliObj = sysCaliDataClass('');
        bg_file_path = [];
        nWaveformToPlot = 100;
        nStep; % plot every nStep measurement

        dataLoadedFlag = 0; % data loaded flag
        rawDataObj % Description
        dataPlotHandle % handle for data lines
        BgPlotHandle % Description
        DataIgnoreLowLineCh1Handle % Description
        DataIgnoreLowLineCh2Handle % Description
        DataIgnoreLowLineCh3Handle % Description

        DataIgnoreHighLineCh1Handle % Description
        DataIgnoreHighLineCh2Handle % Description
        DataIgnoreHighLineCh3Handle % Description

        BgLowLineCh1Handle % Description
        BgLowLineCh2Handle % Description
        BgLowLineCh3Handle % Description

        BgHighLineCh1Handle % Description
        BgHighLineCh2Handle % Description
        BgHighLineCh3Handle % Description

        ChDataArray;
        Ch1Obj; % Description
        Ch2Obj; % Description
        Ch3Obj; % Description
        Ch4Obj; % Description

        ChLaguerreArray;
        Ch1LaguerreObj % Description
        Ch2LaguerreObj
        Ch3LaguerreObj
        Ch4LaguerreObj

        deconChannelSelection% Description
        nPointsPerLine % Description
        nLines
        deconResutlObj % Description
        Ch1PeakLineHandle = []; % Description
        Ch2PeakLineHandle = [];
        Ch3PeakLineHandle = [];
        Ch4PeakLineHandle = [];

        FLImDataObj % FLIm data class object
        Ch1DataObj % channel 1 ChannelDataAPD Obj
        Ch2DataObj % channel 1 ChannelDataAPD Obj
        Ch3DataObj % channel 1 ChannelDataAPD Obj

        Ch1Color = [121,0,141]/255; % channel 1 plot color
        Ch2Color = [0,169,255]/255; % channel 2 plot color
        Ch3Color = [129,255,0]/255; % channel 3 plot color

        fileSelection = 1; % default file selection
        saveFolderName % data folder name
        softwareVersion % Description
    end

    %--------------------------------------------------------------------------

    methods (Access = public)

        function plotData(app,type)
            numOfWFtoPlot = app.nWaveformToPlot;
            plotStep = size(app.Ch1DataObj.preProcessedData,2)/numOfWFtoPlot;
            plotStep = floor(plotStep);
            temp = sortData(app.Ch1DataObj, 'ascend');
            tempMin = min(temp(:));
            if isnan(tempMin)
                tempMin = 0;
            end
            plot(app.UIAxes_ch1DataFig,temp(:,1:plotStep:end))
            app.UIAxes_ch1DataFig.YLim = [tempMin-0.1 2];
            app.UIAxes_ch1DataFig.XLim = [0 size(temp,1)];
            temp = sortData(app.Ch2DataObj, 'ascend');
            if isnan(tempMin)
                tempMin = 0;
            end
            plot(app.UIAxes_ch2DataFig,temp(:,1:plotStep:end))
            app.UIAxes_ch2DataFig.YLim = [tempMin-0.1 2];
            app.UIAxes_ch2DataFig.XLim = [0 size(temp,1)];
            temp = sortData(app.Ch3DataObj, 'ascend');
            if isnan(tempMin)
                tempMin = 0;
            end
            plot(app.UIAxes_ch3DataFig,temp(:,1:plotStep:end))
            app.UIAxes_ch3DataFig.YLim = [tempMin-0.1 2];
            app.UIAxes_ch3DataFig.XLim = [0 size(temp,1)];
            switch type
                case 'raw'
                    titleText1 = ['Ch1 raw FLIm data (no averaging) ' num2str(numOfWFtoPlot) ' waveforms'];
                    titleText2 = ['Ch2 raw FLIm data (no averaging) ' num2str(numOfWFtoPlot) ' waveforms'];
                    titleText3 = ['Ch3 raw FLIm data (no averaging) ' num2str(numOfWFtoPlot) ' waveforms'];
                case 'filtered'
                    titleText1 = ['Ch1 filtered raw FLIm data (no averaging) ' num2str(numOfWFtoPlot) ' waveforms'];
                    titleText2 = ['Ch2 filtered raw FLIm data (no averaging) ' num2str(numOfWFtoPlot) ' waveforms'];
                    titleText3 = ['Ch3 filtered raw FLIm data (no averaging) ' num2str(numOfWFtoPlot) ' waveforms'];
                case 'averaged'
                    titleText1 = ['Ch1 averaged FLIm data (' num2str(app.Ch1DataObj.dataAveraging) ' averaging) ' num2str(numOfWFtoPlot) ' waveforms'];
                    titleText2 = ['Ch2 averaged FLIm data (' num2str(app.Ch2DataObj.dataAveraging) ' averaging) ' num2str(numOfWFtoPlot) ' waveforms'];
                    titleText3 = ['Ch3 averaged FLIm data (' num2str(app.Ch3DataObj.dataAveraging) ' averaging) ' num2str(numOfWFtoPlot) ' waveforms'];
            end
            title(app.UIAxes_ch1DataFig,titleText1)
            title(app.UIAxes_ch2DataFig,titleText2)
            title(app.UIAxes_ch3DataFig,titleText3)

        end


        function plotDataIgnoreMarker(app) % not in use
            if ~isempty(app.DataIgnoreLowLineCh1Handle)
                delete(app.DataIgnoreLowLineCh1Handle)
            end
            if ~isempty(app.DataIgnoreLowLineCh2Handle)
                delete(app.DataIgnoreLowLineCh2Handle)
            end
            if ~isempty(app.DataIgnoreLowLineCh3Handle)
                delete(app.DataIgnoreLowLineCh3Handle)
            end
            if ~isempty(app.DataIgnoreHighLineCh1Handle)
                delete(app.DataIgnoreHighLineCh1Handle)
            end
            if ~isempty(app.DataIgnoreHighLineCh2Handle)
                delete(app.DataIgnoreHighLineCh2Handle)
            end
            if ~isempty(app.DataIgnoreHighLineCh3Handle)
                delete(app.DataIgnoreHighLineCh3Handle)
            end

            lowIdx = app.dataIgnoreLowEditField.Value;
            highIdx = app.dataIgnoreHighEditField.Value;
            app.DataIgnoreLowLineCh1Handle = line(app.UIAxes_ch1BgRemoveFig,[lowIdx lowIdx],[app.UIAxes_ch1BgRemoveFig.YLim],'Color','g','LineStyle','--','LineWidth',1); % plot on data preview figure
            app.DataIgnoreHighLineCh1Handle = line(app.UIAxes_ch1BgRemoveFig,[highIdx highIdx],[app.UIAxes_ch1BgRemoveFig.YLim],'Color','g','LineStyle','--','LineWidth',1) ; %  plot on data preview figure

            app.DataIgnoreLowLineCh2Handle = line(app.UIAxes_ch2BgRemoveFig,[lowIdx lowIdx],[app.UIAxes_ch2BgRemoveFig.YLim],'Color','g','LineStyle','--','LineWidth',1); % plot on BG figure
            app.DataIgnoreHighLineCh2Handle = line(app.UIAxes_ch2BgRemoveFig,[highIdx highIdx],[app.UIAxes_ch2BgRemoveFig.YLim],'Color','g','LineStyle','--','LineWidth',1);  % plot on BG figure

            app.DataIgnoreLowLineCh3Handle = line(app.UIAxes_ch3BgRemoveFig,[lowIdx lowIdx],[app.UIAxes_ch3BgRemoveFig.YLim],'Color','g','LineStyle','--','LineWidth',1); % plot on data preview figure
            app.DataIgnoreHighLineCh3Handle = line(app.UIAxes_ch3BgRemoveFig,[highIdx highIdx],[app.UIAxes_ch3BgRemoveFig.YLim],'Color','g','LineStyle','--','LineWidth',1) ; %  plot on data preview figure
        end

        function plotBGMarker(app)

            if ~isempty(app.BgLowLineCh1Handle)
                delete(app.BgLowLineCh1Handle)
            end
            if ~isempty(app.BgLowLineCh2Handle)
                delete(app.BgLowLineCh2Handle)
            end
            if ~isempty(app.BgLowLineCh3Handle)
                delete(app.BgLowLineCh3Handle)
            end

            if ~isempty(app.BgHighLineCh1Handle)
                delete(app.BgHighLineCh1Handle)
            end
            if ~isempty(app.BgHighLineCh2Handle)
                delete(app.BgHighLineCh2Handle)
            end
            if ~isempty(app.BgHighLineCh3Handle)
                delete(app.BgHighLineCh3Handle)
            end


            lowIdx = app.BgLowEditField.Value;
            highIdx = app.BgHighEditField.Value;
            app.BgLowLineCh1Handle = line(app.UIAxes_ch1DataFig,[lowIdx lowIdx],[app.UIAxes_ch1DataFig.YLim],'Color','r','LineStyle','--','LineWidth',1); % plot on data preview figure
            app.BgHighLineCh1Handle = line(app.UIAxes_ch1DataFig,[highIdx highIdx],[app.UIAxes_ch1DataFig.YLim],'Color','r','LineStyle','--','LineWidth',1) ; %  plot on data preview figure

            app.BgLowLineCh2Handle = line(app.UIAxes_ch2DataFig,[lowIdx lowIdx],[app.UIAxes_ch2DataFig.YLim],'Color','r','LineStyle','--','LineWidth',1); % plot on data preview figure
            app.BgHighLineCh2Handle = line(app.UIAxes_ch2DataFig,[highIdx highIdx],[app.UIAxes_ch2DataFig.YLim],'Color','r','LineStyle','--','LineWidth',1) ; %  plot on data preview figure

            app.BgLowLineCh3Handle = line(app.UIAxes_ch3DataFig,[lowIdx lowIdx],[app.UIAxes_ch3DataFig.YLim],'Color','r','LineStyle','--','LineWidth',1); % plot on data preview figure
            app.BgHighLineCh3Handle = line(app.UIAxes_ch3DataFig,[highIdx highIdx],[app.UIAxes_ch3DataFig.YLim],'Color','r','LineStyle','--','LineWidth',1) ; %  plot on data preview figure

            drawnow
        end

        function plotPreprocessedData(app)
            numOfWFtoPlot = app.nWaveformToPlot;
            ChWidth = app.channelWidthEditField.Value;
            ChWidthInPoints = ChWidth/app.Ch1DataObj.dtUp;
            plotStep = size(app.Ch1DataObj.preProcessedData,2)/numOfWFtoPlot;
            plotStep = floor(plotStep);
            temp = sortData(app.Ch1DataObj, 'ascend');
            tempMin = min(temp(:));
            if isnan(tempMin)
                tempMin=0;
            end
            [~,maxIdx] = max(temp);
            maxIdx = median(maxIdx);
            plot(app.UIAxes_ch1BgRemoveFig,temp(:,1:plotStep:end));
            xline(app.UIAxes_ch1BgRemoveFig,maxIdx-round(0.1*ChWidthInPoints),'c--','LineWidth',2);
            xline(app.UIAxes_ch1BgRemoveFig,maxIdx+round(0.9*ChWidthInPoints),'c--','LineWidth',2);
            yline(app.UIAxes_ch1BgRemoveFig,app.digitizerNoise,'m--','LineWidth',2);
            ylim(app.UIAxes_ch1BgRemoveFig,[tempMin-0.1 2]);
            title(app.UIAxes_ch1BgRemoveFig,'Ch1 DC&BG removed waveforms (100 Waveform)');
            temp = sortData(app.Ch2DataObj, 'ascend');
            tempMin = min(temp(:));
            if isnan(tempMin)
                tempMin=0;
            end
            [~,maxIdx] = max(temp);
            maxIdx = median(maxIdx);
            plot(app.UIAxes_ch2BgRemoveFig,temp(:,1:plotStep:end));
            xline(app.UIAxes_ch2BgRemoveFig,maxIdx-round(0.1*ChWidthInPoints),'c--','LineWidth',2);
            xline(app.UIAxes_ch2BgRemoveFig,maxIdx+round(0.9*ChWidthInPoints),'c--','LineWidth',2);
            yline(app.UIAxes_ch2BgRemoveFig,app.digitizerNoise,'m--','LineWidth',2);
            ylim(app.UIAxes_ch2BgRemoveFig,[tempMin-0.1 2])
            title(app.UIAxes_ch2BgRemoveFig,'Ch2 DC&BG removed waveforms (100 Waveform)')
            temp = sortData(app.Ch3DataObj, 'ascend');
            tempMin = min(temp(:));
            if isnan(tempMin)
                tempMin=0;
            end
            [~,maxIdx] = max(temp);
            maxIdx = median(maxIdx);
            plot(app.UIAxes_ch3BgRemoveFig,temp(:,1:plotStep:end));
            xline(app.UIAxes_ch3BgRemoveFig,maxIdx-round(0.1*ChWidthInPoints),'c--','LineWidth',2);
            xline(app.UIAxes_ch3BgRemoveFig,maxIdx+round(0.9*ChWidthInPoints),'c--','LineWidth',2);
            yline(app.UIAxes_ch3BgRemoveFig,app.digitizerNoise,'m--','LineWidth',2);
            ylim(app.UIAxes_ch3BgRemoveFig,[tempMin-0.1 2]);
            title(app.UIAxes_ch3BgRemoveFig,'Ch3 DC&BG removed waveforms (100 Waveform)')
        end

        function plotUpsampledData(app)
            numOfWFtoPlot = app.nWaveformToPlot;
            plotStep = size(app.Ch1DataObj.preProcessedData,2)/numOfWFtoPlot;
            plotStep = floor(plotStep);
            temp = sortData(app.Ch1DataObj, 'ascend');
            plot(app.UIAxes_ch1BgRemoveFig,temp(:,1:plotStep:end));
            ylim(app.UIAxes_ch1BgRemoveFig,[-0.1 2])
            title(app.UIAxes_ch1BgRemoveFig,'Ch1 DC&BG removed upsampled waveforms (100 Waveform)')
            temp = sortData(app.Ch2DataObj, 'ascend');
            plot(app.UIAxes_ch2BgRemoveFig,temp(:,1:plotStep:end));
            ylim(app.UIAxes_ch2BgRemoveFig,[-0.1 2])
            title(app.UIAxes_ch2BgRemoveFig,'Ch2 DC&BG removed upsampled waveforms (100 Waveform)')
            temp = sortData(app.Ch3DataObj, 'ascend');
            plot(app.UIAxes_ch3BgRemoveFig,temp(:,1:plotStep:end));
            ylim(app.UIAxes_ch3BgRemoveFig,[-0.1 2])
            title(app.UIAxes_ch3BgRemoveFig,'Ch3 DC&BG removed upsampled waveforms (100 Waveform)')
        end

        function plotTruncatedData(app)
            numOfWFtoPlot = app.nWaveformToPlot;
            plotStep = size(app.Ch1DataObj.dataT,2)/numOfWFtoPlot;
            plotStep = floor(plotStep);
            temp = sortTruncatedData(app.Ch1DataObj, 'ascend');
            plot(app.UIAxes_ch1BgRemoveFig,temp(:,1:plotStep:end));
            yline(app.UIAxes_ch1BgRemoveFig,app.digitizerNoise,'m--','LineWidth',2);
            ylim(app.UIAxes_ch1BgRemoveFig,[-0.1 2])
            title(app.UIAxes_ch1BgRemoveFig,'Ch1 DC&BG removed truncated waveforms (100 Waveform)')
            temp = sortTruncatedData(app.Ch2DataObj, 'ascend');
            plot(app.UIAxes_ch2BgRemoveFig,temp(:,1:plotStep:end));
            yline(app.UIAxes_ch2BgRemoveFig,app.digitizerNoise,'m--','LineWidth',2);
            ylim(app.UIAxes_ch2BgRemoveFig,[-0.1 2])
            title(app.UIAxes_ch2BgRemoveFig,'Ch2 DC&BG removed truncated waveforms (100 Waveform)')
            temp = sortTruncatedData(app.Ch3DataObj, 'ascend');
            plot(app.UIAxes_ch3BgRemoveFig,temp(:,1:plotStep:end));
            yline(app.UIAxes_ch3BgRemoveFig,app.digitizerNoise,'m--','LineWidth',2);
            ylim(app.UIAxes_ch3BgRemoveFig,[-0.1 2])
            title(app.UIAxes_ch3BgRemoveFig,'Ch3 DC&BG removed truncated waveforms (100 Waveform)')
        end

        function plotPeakLocations(app,peaks)
            if ~isempty(app.Ch1PeakLineHandle)
                delete (app.Ch1PeakLineHandle)
            end
            if ~isempty(app.Ch2PeakLineHandle)
                delete (app.Ch2PeakLineHandle)
            end
            if ~isempty(app.Ch3PeakLineHandle)
                delete (app.Ch3PeakLineHandle)
            end
            if ~isempty(app.Ch4PeakLineHandle)
                delete (app.Ch4PeakLineHandle)
            end
            app.Ch1PeakLineHandle = line(app.UIAxes_ch1BgRemoveFig,[peaks(1) peaks(1)],[app.UIAxes_ch1DataFig.YLim],'Color','m','LineStyle','--','LineWidth',2);
            app.Ch2PeakLineHandle = line(app.UIAxes_ch1BgRemoveFig,[peaks(2) peaks(2)],[app.UIAxes_ch1DataFig.YLim],'Color','g','LineStyle','--','LineWidth',2);
            app.Ch3PeakLineHandle = line(app.UIAxes_ch1BgRemoveFig,[peaks(3) peaks(3)],[app.UIAxes_ch1DataFig.YLim],'Color','y','LineStyle','--','LineWidth',2);
            app.Ch4PeakLineHandle = line(app.UIAxes_ch1BgRemoveFig,[peaks(4) peaks(4)],[app.UIAxes_ch1DataFig.YLim],'Color','r','LineStyle','--','LineWidth',2);

        end

        function plotChannelData(app,axes,chObj)

            pData = plot(axes,chObj.data(:,1:app.nStep:end));
            hold(axes,'on')
            pIRF = plot(axes,max(max(chObj.data(:,1:app.nStep:end)))/max(chObj.iIRF).*chObj.iIRF,'m','LineWidth',2);
            hold(axes,'off')
            for i = 1:length(pData)
                set(get(get(pData(i),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
            end
            legend(axes,{'IRF'})

        end

        function plotDeconResult(app)
            %get channel number
            channel = app.channelDropDown.Value;
            %get laguerre object
            switch channel
                case 'Channel 1'
                    plotObj = app.Ch1DataObj;
                case 'Channel 2'
                    plotObj = app.Ch2DataObj;
                case 'Channel 3'
                    plotObj = app.Ch3DataObj;
            end

            %get line number
            app.nLines = 2;
            app.nPointsPerLine = size(plotObj.dataT,2);
            % setup points slider bar
            app.LineSlider.Limits = [1 app.nPointsPerLine];
            app.LineSlider.MajorTicks = round(linspace(1, app.nPointsPerLine,10));

            pointIdx = app.LineSlider.Value;

            % calculate linear index
            Idx = pointIdx;
            % get goodness of fit and set line color
            try
                gof = plotObj.stat_test.autoCorr.h(Idx); % goodness of fit to filter data
            catch
                gof=1;
            end
            if gof
                rawColor = [0 0 0 0.5]; % if fit is good, use blue
            else
                rawColor = [0.5 0.5 0.5];% if fit is not good, use grey
            end

            %             autoCorrToPlot = LaguerreObj.stat_test.autoCorr.autoCorrCurve{Idx};
            rawDataToPlot = get(plotObj,'wf_aligned',Idx);
            rawDataMax = max(rawDataToPlot);
            fitToPlotLG = get(plotObj,'fit',Idx);
            residueToPlotLG = rawDataToPlot-fitToPlotLG;
            if app.expDeconFlag
                fitToPlotExp = get(plotObj.expDeconObj,'fit',Idx);
                residueToPlotExp = rawDataToPlot-fitToPlotExp;
            else
                fitToPlotExp = zeros(size(fitToPlotLG));
                residueToPlotExp = zeros(size(residueToPlotLG));
            end

            gain = plotObj.gain(Idx);
            CtrlV = plotObj.CtrlV(Idx);
            LT_LG =  plotObj.Lg_LTs(Idx);
            irfIdx = plotObj.irfIdx(Idx);

            irf = plotObj.APDObj.irfTNorm(:,irfIdx);

            t = plotObj.dtUp*(1:length(rawDataToPlot))';
            %             app.xMaxEditField.Value = max(t); % set maximum x aixs value after loading data
            % plot fitting
            plot(app.UIAxesFitting,t,rawDataToPlot,'.-','Color', rawColor,'LineWidth', 1, 'MarkerSize',10)
            hold(app.UIAxesFitting,'on')
            plot(app.UIAxesFitting,(1:length(irf))*plotObj.dtUp,rawDataMax/max(irf)*irf,'Color',[0 0 1 0.5],'LineWidth', 1)
            plot(app.UIAxesFitting,t,fitToPlotLG,'m','LineWidth', 1.2)
            plot(app.UIAxesFitting,t,fitToPlotExp,'g-.','LineWidth', 1.2)
            hold(app.UIAxesFitting,'off')

            lgd = legend(app.UIAxesFitting,'Raw Data', 'iRF', 'LG fit','mExp fit');
            lgd.FontSize = 8;
            title_temp = sprintf('Point %2d, Laguerre Lifetime %3.2f, Gain %4.0f.',Idx,LT_LG,gain);
            title_temp = [channel title_temp];
            title(app.UIAxesFitting,title_temp)
            ylim(app.UIAxesFitting,[-0.4 plotObj.dataHigh+0.1])
            xlim(app.UIAxesFitting,[app.xMinEditField.Value app.xMaxEditField.Value])


            % plot normalized residue
            normresidueToPlotLG = residueToPlotLG./rawDataMax*100;
            normresidueToPlotExp = residueToPlotExp./rawDataMax*100;
            plot(app.UIAxesResidue,t,normresidueToPlotLG,'r.-', 'LineWidth', 1,'MarkerSize',10);
            hold(app.UIAxesResidue,'on')
            plot(app.UIAxesResidue,t,normresidueToPlotExp,'g.-', 'LineWidth', 1,'MarkerSize',10);
            hold(app.UIAxesResidue,'off')
            lgd = legend(app.UIAxesResidue,'Laguerre', 'Exp','Location','northeast');
            lgd.FontSize = 8;
            try
                ylim(app.UIAxesResidue,[-max(max(norm_residue_to_plot), abs(min(norm_residue_to_plot))) max(max(norm_residue_to_plot), abs(min(norm_residue_to_plot)))]);
            end
            xlim(app.UIAxesResidue,[app.xMinEditField.Value app.xMaxEditField.Value])
            title(app.UIAxesResidue,'Residual (%)')

            %plot auto correlation
            autoCorrToPlotLG = xcorr(residueToPlotLG, 'coeff');
            autoCorrToPlotExp = xcorr(residueToPlotExp, 'coeff');
            plot(app.UIAxesAutoCo, t, autoCorrToPlotLG(length(t):end), 'r.-', 'LineWidth', 2)
            hold(app.UIAxesAutoCo,'on')
            plot(app.UIAxesAutoCo, t, autoCorrToPlotExp(length(t):end), 'g.-', 'LineWidth', 2)
            % setup upper and lower bound of auto-correlation
            bounds(1) = 2 / sqrt(length(rawDataToPlot));
            bounds(2) = -bounds(1);
            plot(app.UIAxesAutoCo, t, (bounds(1) * ones(size(t))), 'k--')
            plot(app.UIAxesAutoCo, t, (bounds(2) * ones(size(t))), 'k--')
            hold(app.UIAxesAutoCo,'off')
            ylabel(app.UIAxesAutoCo,'Res. Autocorrelation (a.u.)');
            try
                ylim([-max(max(autocorr_to_plot), abs(min(autocorr_to_plot))) max(max(autocorr_to_plot), abs(min(autocorr_to_plot)))]);
            end
            xlim(app.UIAxesAutoCo,[app.xMinEditField.Value app.xMaxEditField.Value])
            lgd = legend(app.UIAxesAutoCo,'Laguerre', 'Exp','Location','northeast');
            lgd.FontSize = 8;
        end

        function plotImg(app)
            %get channel number
            channel = app.deconChannelSelection(str2double(app.channelDropDown.Value));
            %get laguerre object

            imgLT = app.deconResutlObj.ltMap{channel};
            imgLT = imgLT;
            imagesc(app.UIAxesImg,imgLT)
            hold(app.UIAxesImg,'on')
            plot(app.UIAxesImg,app.LineSliderSlider.Value,app.LineSlider.Value,'w+', 'MarkerSize', 20,'LineWidth',3);
            hold(app.UIAxesImg,'off')
            xlim(app.UIAxesImg,[1 app.nPointsPerLine])
            ylim(app.UIAxesImg,[1 app.nLines])
            %             axis(app.UIAxesImg,'equal')
            colormap(app.UIAxesImg,'jet')
            colorbar(app.UIAxesImg)

        end


        function updateChannelDropDown(app)
            dropDownList = cell(1);
            if app.channel1CheckBox.Value
                dropDownList = [dropDownList; app.channel1CheckBox.Text];
            end
            if app.channel2CheckBox.Value
                dropDownList = [dropDownList; app.channel2CheckBox.Text];
            end
            if app.channel3CheckBox.Value
                dropDownList = [dropDownList; app.channel3CheckBox.Text];
            end
            app.channelDropDown.Items = dropDownList(2:end);
        end

        function clearAllPreviewAxis(app)
            cla(app.UIAxes_ch1DataFig)
            cla(app.UIAxes_ch2DataFig)
            cla(app.UIAxes_ch3DataFig)
            cla(app.UIAxes_ch1BgRemoveFig)
            cla(app.UIAxes_ch2BgRemoveFig)
            cla(app.UIAxes_ch3BgRemoveFig)
        end

        function saveAPDDetectorPath(app)
            ini = IniConfig();
            ini.AddSections({'APD Detector Temp Path'});
            ini.AddKeys('APD Detector Temp Path', 'ADP1', app.apd1FileEditField.Value);
            ini.AddKeys('APD Detector Temp Path', 'ADP2', app.apd2FileEditField.Value);
            ini.AddKeys('APD Detector Temp Path', 'ADP3', app.apd3FileEditField.Value);
            ini.WriteFile(fullfile(app.rootFolder,'Algorithms','Preprocessing_APD','DetectorPathTamp.ini'));
        end

        function saveDeconSettingIni(app)
            ini = IniConfig();
            ini.AddSections({'DeCon Parameter'});
            ini.AddKeys('DeCon Parameter', 'Time Resolution (ns)', app.TimeResolutionnsEditField.Value);
            ini.AddKeys('DeCon Parameter', 'Channel Width (ns)', app.channelWidthEditField.Value);
            ini.AddKeys('DeCon Parameter', 'Signal Low (V)', app.LowThresholdEditField.Value);
            ini.AddKeys('DeCon Parameter', 'Signal High (V)', app.HighThresholdEditField.Value);
            ini.AddKeys('DeCon Parameter', 'Laguerre Order', app.LaguerreOrderEditField.Value);
            ini.AddKeys('DeCon Parameter', 'Ch1 Alpha Value', app.Ch1AlphaValueDropDown.Value);
            ini.AddKeys('DeCon Parameter', 'Ch2 Alpha Value', app.Ch2AlphaValueDropDown.Value);
            ini.AddKeys('DeCon Parameter', 'Ch3 Alpha Value', app.Ch3AlphaValueDropDown.Value);
            ini.AddKeys('DeCon Parameter', 'Deconvolve Channel 1?', int64(app.channel1CheckBox.Value));
            ini.AddKeys('DeCon Parameter', 'Deconvolve Channel 2?', int64(app.channel2CheckBox.Value));
            ini.AddKeys('DeCon Parameter', 'Deconvolve Channel 3?', int64(app.channel3CheckBox.Value));
            ini.AddKeys('DeCon Parameter', 'BG Index Low', app.BgLowEditField.Value);
            ini.AddKeys('DeCon Parameter', 'BG Index High', app.BgHighEditField.Value);
            ini.AddKeys('DeCon Parameter', 'Data Ignore Index Low', app.dataIgnoreLowEditField.Value);
            ini.AddKeys('DeCon Parameter', 'Data Ignore Index High', app.dataIgnoreHighEditField.Value);
            ini.WriteFile(fullfile(app.rootFolder,'Algorithms','Preprocessing_APD','DeConSetting.ini'));
        end

        function setAPDDetectorPath(app)
            %-------------------read APD detector file---------------------
            ini = IniConfig();
            ini.ReadFile('DetectorPathTamp.ini');
            sections = ini.GetSections();
            [keys, ~] = ini.GetKeys(sections{1});
            APDPath = ini.GetValues(sections{1}, keys);
            % load ch1 APD file
            try
                app.apd1Obj = apdClass(APDPath{1});
                creatFromPath(app.apd1Obj);
                app.apd1LoadedLamp.Color = 'g';
                app.apd1LoadedFlag = 1;
                app.apd1FileEditField.Value = app.apd1Obj.filePath;
                app.apd1ViewButton.Enable = 'On'; % enable view button
                [filepath,name,ext] = fileparts(APDPath{1});
                app.dataInfoObj.apd1Name = [name ext];
                app.dataInfoObj.apd1Folder = filepath;
            catch
                warndlg('Could not load channel 1 APD detector file, please manually by clicking load button');
                app.apd1FileEditField.Value = 'Please Select Channel 1 APD File!';
            end
            %load ch2 APD file
            try
                app.apd2Obj = apdClass(APDPath{2});
                creatFromPath(app.apd2Obj); % load 1st, if error throw warning
                app.apd2LoadedLamp.Color = 'g';
                app.apd2LoadedFlag = 1;
                app.apd2FileEditField.Value = app.apd2Obj.filePath;
                app.apd2ViewButton.Enable = 'On'; % enable view button
                [filepath,name,ext] = fileparts(APDPath{2}); %
                app.dataInfoObj.apd2Name = [name ext];
                app.dataInfoObj.apd2Folder = filepath;
            catch
                warndlg('Could not load channel 2 APD detector file, please manually by clicking load button');
                app.apd2FileEditField.Value = 'Please Select Channel 2 APD File!';
            end

            try
                app.apd3Obj = apdClass(APDPath{3});
                creatFromPath(app.apd3Obj);
                app.apd3LoadedLamp.Color = 'g';
                app.apd3LoadedFlag = 1;
                app.apd3FileEditField.Value = app.apd3Obj.filePath;
                app.apd3ViewButton.Enable = 'On'; % enable view button
                [filepath,name,ext] = fileparts(APDPath{3}); %
                app.dataInfoObj.apd3Name = [name ext];
                app.dataInfoObj.apd3Folder = filepath;
            catch
                warndlg('Could not load channel 3 APD detector file, please manually by clicking load button');
                app.apd3FileEditField.Value = 'Please Select Channel 3 APD File!';
            end
            %--------------------------------------------------------------
        end

        function loadDeConSettingIni(app)
            ini = IniConfig();
            ini.ReadFile('DeConSetting.ini');
            sections = ini.GetSections();
            [keys, ~] = ini.GetKeys(sections{1});
            DeConSettings = ini.GetValues(sections{1}, keys);
            app.TimeResolutionnsEditField.Value = DeConSettings{1};
            app.channelWidthEditField.Value = DeConSettings{2};
            app.LowThresholdEditField.Value = DeConSettings{3};
            app.HighThresholdEditField.Value = DeConSettings{4};
            app.LaguerreOrderEditField.Value = DeConSettings{5};
            app.Ch1AlphaValueDropDown.Value = num2str(DeConSettings{6},'%.3f');
            app.Ch2AlphaValueDropDown.Value = num2str(DeConSettings{7},'%.3f');
            app.Ch3AlphaValueDropDown.Value = num2str(DeConSettings{8},'%.3f');
            app.channel1CheckBox.Value = DeConSettings{9};
            app.channel2CheckBox.Value = DeConSettings{10};
            app.channel3CheckBox.Value = DeConSettings{11};
            app.BgLowEditField.Value = DeConSettings{12};
            app.BgHighEditField.Value = DeConSettings{13};
            app.dataIgnoreLowEditField.Value = DeConSettings{14};
            app.dataIgnoreHighEditField.Value = DeConSettings{15};
        end

        function save_LG_vs_mExp_lifetime_fig(app)
            channel = app.channelDropDown.Value;
            switch channel
                case 'Channel 1'
                    plotObj = app.Ch1DataObj;
                    %                             c = app.Ch1Color;
                    tl = 'Channel 1';
                case 'Channel 2'
                    plotObj = app.Ch2DataObj;
                    %                             c = app.Ch2Color;
                    tl = 'Channel 2';
                case 'Channel 3'
                    plotObj = app.Ch3DataObj;
                    %                             c = app.Ch3Color;
                    tl = 'Channel 3';
            end
            lg = {'LG decay vs mExp decay','LG decay vs mExp formula','y=x'};
            LT_LG_decay = plotObj.Lg_LTs;
            LT_mExp_decay = plotObj.expDeconObj.LTs_decay;
            LT_mExp_formula = plotObj.expDeconObj.LTs_formula;
            f1 = figure('Position',[100 100 960 540]);
            ax = gca;
            scatter(ax,LT_LG_decay,LT_mExp_decay, [],'green','filled','o');
            hold(ax,'on');
            scatter(ax,LT_LG_decay,LT_mExp_formula, [],'m','filled','s');
            plot(ax,[0 20],[0 20],'r--','LineWidth',2)
            hold(ax,'off');
            %                     axis(app.UIAxesDeconResult,'tight')
            xlabel(ax,'Laguerre lifetimes (ns)')
            ylabel(ax,'Muti-exponential Lifetimes (ns)')
            title(ax,[tl ' Laguerre vs Multi-exponential lifetimes'])
            legend(ax,lg,'Location','northwest')
            xlim(ax,[mean(LT_LG_decay,'omitnan')-3*std(LT_LG_decay,'omitnan') mean(LT_LG_decay,'omitnan')+3*std(LT_LG_decay,'omitnan')])
            ylim(ax,[mean([LT_mExp_decay;LT_mExp_decay],'omitnan')-3*std([LT_mExp_decay;LT_mExp_decay],'omitnan') mean([LT_mExp_decay;LT_mExp_decay],'omitnan')+3*std([LT_mExp_decay;LT_mExp_decay],'omitnan')])
            filename = fullfile(app.dataFolder, app.saveFolderName,'LG_vs_mEXp_Lifetime.fig');
            saveas(f1,filename)
            close(f1)
        end
    end

    methods (Access = private)

        function uiUpdate(app,step)
            switch step
                case 'init'
                    % diable
                    app.LowThresholdEditField.Enable = 0;
                    app.HighThresholdEditField.Enable = 0;
                    app.DataAverageDropDown.Enable = 0;
                    app.BgThreholdEditField.Enable = 0;
                    app.BgHighEditField.Enable = 0;
                    app.BgLowEditField.Enable = 0;
                    app.dataIgnoreLowEditField.Enable = 0;
                    app.dataIgnoreHighEditField.Enable = 0;
                    app.TimeResolutionnsEditField.Enable = 0;
                    app.channelWidthEditField.Enable = 0;
                    app.channel1CheckBox.Enable = 0;
                    app.channel2CheckBox.Enable = 0;
                    app.channel3CheckBox.Enable = 0;
                    app.LaguerreOrderEditField.Enable = 0;
                    app.Ch1AlphaValueDropDown.Enable = 0;
                    app.Ch2AlphaValueDropDown.Enable = 0;
                    app.Ch3AlphaValueDropDown.Enable = 0;
                    app.RemoveSaturationButton.Enable = 0;
                    app.AverageDataButton.Enable = 0;
                    app.RemoveBGButton.Enable = 0;
                    app.UpsampleButton.Enable = 0;
                    app.TruncateButton.Enable = 0;
                    app.DeconButton.Enable = 0;
                    app.RunEXPCheckBox.Enable = 0;
                    app.exponentialsDropDown.Enable = 0;
                    % enable
                    app.LoadDataButton.Enable = 1;

                case 'dataLoaded'
                    % diable
                    app.DataAverageDropDown.Enable = 0;
                    app.BgThreholdEditField.Enable = 0;
                    app.BgHighEditField.Enable = 0;
                    app.BgLowEditField.Enable = 0;
                    app.dataIgnoreLowEditField.Enable = 0;
                    app.dataIgnoreHighEditField.Enable = 0;
                    app.TimeResolutionnsEditField.Enable = 0;
                    app.channelWidthEditField.Enable = 0;
                    app.channel1CheckBox.Enable = 0;
                    app.channel2CheckBox.Enable = 0;
                    app.channel3CheckBox.Enable = 0;
                    app.LaguerreOrderEditField.Enable = 0;
                    app.Ch1AlphaValueDropDown.Enable = 0;
                    app.Ch2AlphaValueDropDown.Enable = 0;
                    app.Ch3AlphaValueDropDown.Enable = 0;
                    app.AverageDataButton.Enable = 0;
                    app.RemoveBGButton.Enable = 0;
                    app.TruncateButton.Enable = 0;
                    app.DeconButton.Enable = 0;
                    app.LowThresholdEditField.Enable = 0;
                    app.HighThresholdEditField.Enable = 0;
                    app.RemoveSaturationButton.Enable = 0;
                    % enable
                    app.UpsampleButton.Enable = 1;
                    app.LoadDataButton.Enable = 1;

                case 'Upsampled'
                    % diable
                    app.LowThresholdEditField.Enable = 0;
                    app.HighThresholdEditField.Enable = 0;
                    app.BgThreholdEditField.Enable = 0;
                    app.BgHighEditField.Enable = 0;
                    app.BgLowEditField.Enable = 0;
                    app.dataIgnoreLowEditField.Enable = 0;
                    app.dataIgnoreHighEditField.Enable = 0;
                    app.TimeResolutionnsEditField.Enable = 0;
                    app.channelWidthEditField.Enable = 0;
                    app.channel1CheckBox.Enable = 0;
                    app.channel2CheckBox.Enable = 0;
                    app.channel3CheckBox.Enable = 0;
                    app.LaguerreOrderEditField.Enable = 0;
                    app.Ch1AlphaValueDropDown.Enable = 0;
                    app.Ch2AlphaValueDropDown.Enable = 0;
                    app.Ch3AlphaValueDropDown.Enable = 0;
                    app.AverageDataButton.Enable = 0;
                    app.RemoveBGButton.Enable = 0;
                    app.UpsampleButton.Enable = 0;
                    app.TruncateButton.Enable = 0;
                    app.DeconButton.Enable = 0;
                    app.DataAverageDropDown.Enable = 0;
                    % enable
                    app.LoadDataButton.Enable = 1;
                    app.RemoveSaturationButton.Enable = 1;

                case 'removedSat'
                    % diable
                    app.LowThresholdEditField.Enable = 0;
                    app.HighThresholdEditField.Enable = 0;
                    app.RemoveSaturationButton.Enable = 0;
                    app.BgThreholdEditField.Enable = 0;
                    app.BgHighEditField.Enable = 0;
                    app.BgLowEditField.Enable = 0;
                    app.dataIgnoreLowEditField.Enable = 0;
                    app.dataIgnoreHighEditField.Enable = 0;
                    app.TimeResolutionnsEditField.Enable = 0;
                    app.channelWidthEditField.Enable = 0;
                    app.channel1CheckBox.Enable = 0;
                    app.channel2CheckBox.Enable = 0;
                    app.channel3CheckBox.Enable = 0;
                    app.LaguerreOrderEditField.Enable = 0;
                    app.Ch1AlphaValueDropDown.Enable = 0;
                    app.Ch2AlphaValueDropDown.Enable = 0;
                    app.Ch3AlphaValueDropDown.Enable = 0;
                    app.RemoveBGButton.Enable = 0;
                    app.UpsampleButton.Enable = 0;
                    app.TruncateButton.Enable = 0;
                    app.DeconButton.Enable = 0;
                    % enable
                    app.LoadDataButton.Enable = 1;
                    app.DataAverageDropDown.Enable = 1;
                    app.AverageDataButton.Enable = 1;

                case 'dataAveraged'
                    % diable
                    app.LowThresholdEditField.Enable = 0;
                    app.HighThresholdEditField.Enable = 0;
                    app.RemoveSaturationButton.Enable = 0;
                    app.BgThreholdEditField.Enable = 0;
                    app.TimeResolutionnsEditField.Enable = 0;
                    app.channelWidthEditField.Enable = 0;
                    app.channel1CheckBox.Enable = 0;
                    app.channel2CheckBox.Enable = 0;
                    app.channel3CheckBox.Enable = 0;
                    app.LaguerreOrderEditField.Enable = 0;
                    app.Ch1AlphaValueDropDown.Enable = 0;
                    app.Ch2AlphaValueDropDown.Enable = 0;
                    app.Ch3AlphaValueDropDown.Enable = 0;
                    app.AverageDataButton.Enable = 0;
                    app.UpsampleButton.Enable = 0;
                    app.TruncateButton.Enable = 0;
                    app.DeconButton.Enable = 0;
                    app.DataAverageDropDown.Enable = 0;
                    app.RemoveBGButton.Enable = 1;
                    app.BgHighEditField.Enable = 1;
                    app.BgLowEditField.Enable = 1;
                    app.dataIgnoreLowEditField.Enable = 0;
                    app.dataIgnoreHighEditField.Enable = 0;
                    app.LoadDataButton.Enable = 1;

                case 'BGRemoved'
                    % diable
                    app.LowThresholdEditField.Enable = 0;
                    app.HighThresholdEditField.Enable = 0;
                    app.RemoveSaturationButton.Enable = 0;
                    app.BgThreholdEditField.Enable = 0;
                    app.BgHighEditField.Enable = 0;
                    app.BgLowEditField.Enable = 0;
                    app.dataIgnoreLowEditField.Enable = 0;
                    app.dataIgnoreHighEditField.Enable = 0;
                    app.TimeResolutionnsEditField.Enable = 0;
                    app.channelWidthEditField.Enable = 1;
                    app.channel1CheckBox.Enable = 0;
                    app.channel2CheckBox.Enable = 0;
                    app.channel3CheckBox.Enable = 0;
                    app.LaguerreOrderEditField.Enable = 0;
                    app.Ch1AlphaValueDropDown.Enable = 0;
                    app.Ch2AlphaValueDropDown.Enable = 0;
                    app.Ch3AlphaValueDropDown.Enable = 0;
                    app.AverageDataButton.Enable = 0;
                    app.RemoveBGButton.Enable = 0;
                    app.UpsampleButton.Enable = 0;
                    app.TruncateButton.Enable = 1;
                    app.DeconButton.Enable = 0;
                    app.DataAverageDropDown.Enable = 0;
                    app.exponentialsDropDown.Enable = 0;
                    app.LoadDataButton.Enable = 1;

                case 'truncated'
                    % diable
                    app.LowThresholdEditField.Enable = 0;
                    app.HighThresholdEditField.Enable = 0;
                    app.RemoveSaturationButton.Enable = 0;
                    app.BgThreholdEditField.Enable = 0;
                    app.BgHighEditField.Enable = 0;
                    app.BgLowEditField.Enable = 0;
                    app.dataIgnoreLowEditField.Enable = 1;
                    app.dataIgnoreHighEditField.Enable = 1;
                    app.TimeResolutionnsEditField.Enable = 0;
                    app.channelWidthEditField.Enable = 0;
                    app.channel1CheckBox.Enable = 0;
                    app.channel2CheckBox.Enable = 0;
                    app.channel3CheckBox.Enable = 0;
                    app.LaguerreOrderEditField.Enable = 1;
                    app.Ch1AlphaValueDropDown.Enable = 1;
                    app.Ch2AlphaValueDropDown.Enable = 1;
                    app.Ch3AlphaValueDropDown.Enable = 1;
                    app.AverageDataButton.Enable = 0;
                    app.RemoveBGButton.Enable = 0;
                    app.TruncateButton.Enable = 0;
                    app.DeconButton.Enable = 1;
                    app.DataAverageDropDown.Enable = 0;
                    app.RunEXPCheckBox.Enable = 1;
                    app.exponentialsDropDown.Enable = 0;
                    app.LoadDataButton.Enable = 1;

            end

        end

        function clearAllFittingResultPlot(app)
            cla(app.UIAxesFitting)
            cla(app.UIAxesResidue)
            cla(app.UIAxesAutoCo)
            cla(app.UIAxesDeconResult)
        end

        function plotTrace(app,option)
            cla(app.UIAxesDeconResult)
            switch option
                case 'Decon Intensity'
                    plot(app.UIAxesDeconResult,app.Ch1DataObj.Lg_INTsGainCorrected,'.-', 'Color',app.Ch1Color)
                    hold(app.UIAxesDeconResult,'on')
                    plot(app.UIAxesDeconResult,app.Ch2DataObj.Lg_INTsGainCorrected,'.-', 'Color',app.Ch2Color)
                    plot(app.UIAxesDeconResult,app.Ch3DataObj.Lg_INTsGainCorrected,'.-', 'Color',app.Ch3Color)
                    hold(app.UIAxesDeconResult,'off')
                    legend(app.UIAxesDeconResult,'Ch1','Ch2','Ch3')
                    xlabel(app.UIAxesDeconResult,'Points')
                    ylabel(app.UIAxesDeconResult,'Decon Intensity (a.u.)')
                    title(app.UIAxesDeconResult,'Decon Intensity Trace Decon Data Points Gain Corrected')
                    axis(app.UIAxesDeconResult,'tight')
                case 'Raw Intensity'
                    INT1 = sum(app.Ch1DataObj.preProcessedData)';
                    INT2= sum(app.Ch2DataObj.preProcessedData)';
                    INT3 = sum(app.Ch3DataObj.preProcessedData)';
                    INT1GC = INT1./app.Ch1DataObj.gain;
                    INT2GC = INT2./app.Ch2DataObj.gain;
                    INT3GC = INT3./app.Ch3DataObj.gain;
                    plot(app.UIAxesDeconResult,INT1GC,'.-', 'Color',app.Ch1Color)
                    hold(app.UIAxesDeconResult,'on')
                    plot(app.UIAxesDeconResult,INT2GC,'.-', 'Color',app.Ch2Color)
                    plot(app.UIAxesDeconResult,INT3GC,'.-', 'Color',app.Ch3Color)
                    hold(app.UIAxesDeconResult,'off')
                    axis(app.UIAxesDeconResult,'tight')
                    ylim(app.UIAxesDeconResult,[quantile([INT1GC;INT2GC;INT3GC],0.01) quantile([INT1GC;INT2GC;INT3GC],0.95)])
                    legend(app.UIAxesDeconResult,'Ch1','Ch2','Ch3')
                    xlabel(app.UIAxesDeconResult,'Points')
                    ylabel(app.UIAxesDeconResult,'Raw Intensity (a.u.)')
                    title(app.UIAxesDeconResult,'Raw Intensity Trace Decon Data Points Gain Corrected')
                case 'Gain'
                    G1 = app.Ch1DataObj.gain;
                    G2 = app.Ch2DataObj.gain;
                    G3 = app.Ch3DataObj.gain;
                    plot(app.UIAxesDeconResult,G1,'.-', 'Color',app.Ch1Color)
                    hold(app.UIAxesDeconResult,'on')
                    plot(app.UIAxesDeconResult,G2,'.-', 'Color',app.Ch2Color)
                    plot(app.UIAxesDeconResult,G3,'.-', 'Color',app.Ch3Color)
                    hold(app.UIAxesDeconResult,'off')
                    axis(app.UIAxesDeconResult,'tight')
                    legend(app.UIAxesDeconResult,'Ch1','Ch2','Ch3')
                    xlabel(app.UIAxesDeconResult,'Points')
                    ylabel(app.UIAxesDeconResult,'Gain')
                    title(app.UIAxesDeconResult,'Gain Trace Decon Data Points Gain Corrected')
                case 'Laguerre basis'
                    channel = app.channelDropDown.Value;
                    switch channel
                        case 'Channel 1'
                            plotObj = app.Ch1DataObj;
                        case 'Channel 2'
                            plotObj = app.Ch2DataObj;
                        case 'Channel 3'
                            plotObj = app.Ch3DataObj;
                    end
                    plot(app.UIAxesDeconResult,plotObj.LaguerreBasis)
                    axis(app.UIAxesDeconResult,'tight')
                    title(app.UIAxesDeconResult,'Laguerre Basis')

                case 'Raw Waveforms'
                    channel = app.channelDropDown.Value;
                    switch channel
                        case 'Channel 1'
                            plotObj = app.Ch1DataObj;
                        case 'Channel 2'
                            plotObj = app.Ch2DataObj;
                        case 'Channel 3'
                            plotObj = app.Ch3DataObj;
                    end
                    pointIdx = app.LineSlider.Value;
                    rawDataIdx = (pointIdx-1)*plotObj.dataAveraging+1:pointIdx*plotObj.dataAveraging;
                    rawDatatoPlot = plotObj.rawDataUpsampled(:,rawDataIdx);
                    plot(app.UIAxesDeconResult,rawDatatoPlot,'Marker',"*")
                    axis(app.UIAxesDeconResult,'tight')
                    xlabel(app.UIAxesDeconResult, 'Points')
                    ylabel(app.UIAxesDeconResult, 'Voltage (V)')
                    lg = cell(plotObj.dataAveraging,1);
                    for i = 1:plotObj.dataAveraging
                        temp = sprintf('WF #%d',i);
                        lg{i} = temp;
                    end
                    legend(app.UIAxesDeconResult,lg)
                case 'Residual Max'
                    channel = app.channelDropDown.Value;
                    switch channel
                        case 'Channel 1'
                            plotObj = app.Ch1DataObj;
                            c = app.Ch1Color;
                            lg = 'Channel 1 Residual';
                        case 'Channel 2'
                            plotObj = app.Ch2DataObj;
                            c = app.Ch2Color;
                            lg = 'Channel 2 Residual';
                        case 'Channel 3'
                            plotObj = app.Ch3DataObj;
                            c = app.Ch3Color;
                            lg = 'Channel 3 Residual';
                    end
                    res = get(plotObj,'res')*100;
                    resMax = max(res);
                    plot(app.UIAxesDeconResult, resMax, 'Color', c, 'Marker','*','LineStyle',"none");
                    axis(app.UIAxesDeconResult,'tight')
                    ylim(app.UIAxesDeconResult,[quantile(resMax,0.1) quantile(resMax,0.9)])
                    xlabel(app.UIAxesDeconResult,'Point #')
                    ylabel(app.UIAxesDeconResult,'Max residual (%)')
                    title(app.UIAxesDeconResult,'Max residual vs point measurement #')
                    legend(app.UIAxesDeconResult,lg)
                case 'Lifetime'
                    channel = app.channelDropDown.Value;
                    switch channel
                        case 'Channel 1'
                            plotObj = app.Ch1DataObj;
                            c = app.Ch1Color;
                            lg = 'Channel 1 lifetime';
                        case 'Channel 2'
                            plotObj = app.Ch2DataObj;
                            c = app.Ch2Color;
                            lg = 'Channel 2 lifetime';
                        case 'Channel 3'
                            plotObj = app.Ch3DataObj;
                            c = app.Ch3Color;
                            lg = 'Channel 3 lifetime';
                    end
                    LT = plotObj.Lg_LTs;
                    plot(app.UIAxesDeconResult, LT, 'Marker',"*",'Color',c);
                    axis(app.UIAxesDeconResult,'tight')
                    xlabel(app.UIAxesDeconResult,'Point #')
                    ylabel(app.UIAxesDeconResult,'Lifetime (ns)')
                    title(app.UIAxesDeconResult,'Lifetime vs point measurement #')
                    ylim(app.UIAxesDeconResult,[quantile(LT,0.1) quantile(LT,0.9)])
                    legend(app.UIAxesDeconResult,lg)
                case 'SNR'
                    channel = app.channelDropDown.Value;
                    switch channel
                        case 'Channel 1'
                            plotObj = app.Ch1DataObj;
                            c = app.Ch1Color;
                            lg = 'Channel 1 SNR';
                        case 'Channel 2'
                            plotObj = app.Ch2DataObj;
                            c = app.Ch2Color;
                            lg = 'Channel 2 SNR';
                        case 'Channel 3'
                            plotObj = app.Ch3DataObj;
                            c = app.Ch3Color;
                            lg = 'Channel 3 SNR';
                    end
                    SNR = plotObj.SNR;
                    plot(app.UIAxesDeconResult, SNR, 'Marker',"*",'Color',c);
                    %                     axis(app.UIAxesDeconResult,'tight')
                    xlabel(app.UIAxesDeconResult,'Point #')
                    ylabel(app.UIAxesDeconResult,'SNR (dB)')
                    title(app.UIAxesDeconResult,'SNR vs point measurement #')
                    app.UIAxesDeconResult.YLim = [mean(SNR,'omitnan')-2*std(SNR,'omitnan') mean(SNR,'omitnan')+2*std(SNR,'omitnan')];
                    legend(app.UIAxesDeconResult,lg)
                case 'LG vs mExp'
                    channel = app.channelDropDown.Value;
                    switch channel
                        case 'Channel 1'
                            plotObj = app.Ch1DataObj;
                            %                             c = app.Ch1Color;
                            tl = 'Channel 1';
                        case 'Channel 2'
                            plotObj = app.Ch2DataObj;
                            %                             c = app.Ch2Color;
                            tl = 'Channel 2';
                        case 'Channel 3'
                            plotObj = app.Ch3DataObj;
                            %                             c = app.Ch3Color;
                            tl = 'Channel 3';
                    end
                    lg = {'LG decay vs mExp decay','LG decay vs mExp formula','y=x'};
                    LT_LG_decay = plotObj.Lg_LTs;
                    LT_mExp_decay = plotObj.expDeconObj.LTs_decay;
                    LT_mExp_formula = plotObj.expDeconObj.LTs_formula;
                    axes(app.UIAxesDeconResult)
                    scatter(app.UIAxesDeconResult,LT_LG_decay,LT_mExp_decay, [],'green','filled','o');
                    hold(app.UIAxesDeconResult,'on');
                    scatter(app.UIAxesDeconResult,LT_LG_decay,LT_mExp_formula, [],'m','filled','s');
                    plot(app.UIAxesDeconResult,[0 20],[0 20],'r--','LineWidth',2)
                    hold(app.UIAxesDeconResult,'off');
                    %                     axis(app.UIAxesDeconResult,'tight')
                    xlabel(app.UIAxesDeconResult,'Laguerre lifetimes (ns)')
                    ylabel(app.UIAxesDeconResult,'Muti-exponential Lifetimes (ns)')
                    title(app.UIAxesDeconResult,[tl ' Laguerre vs Multi-exponential lifetimes'])
                    legend(app.UIAxesDeconResult,lg,'Location','northwest')
                    xlim(app.UIAxesDeconResult,[mean(LT_LG_decay,'omitnan')-3*std(LT_LG_decay,'omitnan') mean(LT_LG_decay,'omitnan')+3*std(LT_LG_decay,'omitnan')])
                    ylim(app.UIAxesDeconResult,[mean([LT_mExp_decay;LT_mExp_decay],'omitnan')-3*std([LT_mExp_decay;LT_mExp_decay],'omitnan') mean([LT_mExp_decay;LT_mExp_decay],'omitnan')+3*std([LT_mExp_decay;LT_mExp_decay],'omitnan')])

                case 'Square Error'
                    channel = app.channelDropDown.Value;
                    switch channel
                        case 'Channel 1'
                            plotObj = app.Ch1DataObj;
                            %                             c = app.Ch1Color;
                            tl = 'Channel 1';
                        case 'Channel 2'
                            plotObj = app.Ch2DataObj;
                            %                             c = app.Ch2Color;
                            tl = 'Channel 2';
                        case 'Channel 3'
                            plotObj = app.Ch3DataObj;
                            %                             c = app.Ch3Color;
                            tl = 'Channel 3';
                    end
                    res_Lg = get(plotObj,'wf_aligned')-get(plotObj,'fit');
                    se_Lg = sum(res_Lg.^2)';
                    res_Exp = plotObj.expDeconObj.spec_aligned-get(plotObj.expDeconObj,'fit');
                    se_Exp = sum(res_Exp.^2)';
                    scatter(app.UIAxesDeconResult,se_Lg,se_Exp, [],'magenta','filled','o');
                    hold(app.UIAxesDeconResult,'on');
                    plot(app.UIAxesDeconResult,[0 max(se_Exp)],[0 max(se_Exp)],'g--','LineWidth',1.5)
                    hold(app.UIAxesDeconResult,'off');
                    xlabel(app.UIAxesDeconResult,'Laguerre Square Error')
                    ylabel(app.UIAxesDeconResult,'Muti-exponential Square Error')
                    axis(app.UIAxesDeconResult,'tight')
                    title(app.UIAxesDeconResult,'Square Error')
                    lg = {'Laguerre vs Exponential','y=x'};
                    legend(app.UIAxesDeconResult,lg,'Location','northwest')
                    ylim(app.UIAxesDeconResult,[quantile(se_Exp,0.1) quantile(se_Exp,0.9)])
                    xlim(app.UIAxesDeconResult,[quantile(se_Lg,0.1) quantile(se_Lg,0.9)])

                case 'Shift'
                    channel = app.channelDropDown.Value;
                    switch channel
                        case 'Channel 1'
                            plotObj = app.Ch1DataObj;
                            c = app.Ch1Color;
                            tl = 'Channel 1';
                        case 'Channel 2'
                            plotObj = app.Ch2DataObj;
                            c = app.Ch2Color;
                            tl = 'Channel 2';
                        case 'Channel 3'
                            plotObj = app.Ch3DataObj;
                            c = app.Ch3Color;
                            tl = 'Channel 3';
                    end
                    shift = plotObj.shift;
                    plot(app.UIAxesDeconResult, shift, 'Marker',"*",'Color',c);
                    ylim(app.UIAxesDeconResult,[min(shift)-5 max(shift)+5])
                    %                     axis(app.UIAxesDeconResult,'tight')
                    xlabel(app.UIAxesDeconResult,'Point #')
                    ylabel(app.UIAxesDeconResult,'iRF shift (points)')
                    title(app.UIAxesDeconResult,'optimal iRF shift')

            end
        end


        function clearBGAxis(app)
            cla(app.bg1UIAxes)
            cla(app.bg2UIAxes)
            cla(app.bg3UIAxes)
        end

        function plotBG(app)
            plot(app.bg1UIAxes,app.bgObj.bgCh1,'Marker','.','LineStyle','--','Color', app.Ch1Color,'MarkerSize',8);
            title(app.bg1UIAxes,sprintf('Channel 1 Background, CtrlV = %.3f, Gain = %.3f',app.bgObj.CtrlV1, app.bgObj.bgGain1));
            app.bg1UIAxes.YLim = [min(app.bgObj.bgCh1)-0.1 max(app.bgObj.bgCh1)+0.1];
            plot(app.bg2UIAxes,app.bgObj.bgCh2,'Marker','.','LineStyle','--','Color', app.Ch2Color,'MarkerSize',8);
            title(app.bg2UIAxes,sprintf('Channel 2 Background, CtrlV = %.3f, Gain = %.3f',app.bgObj.CtrlV2, app.bgObj.bgGain2));
            app.bg2UIAxes.YLim = [min(app.bgObj.bgCh2)-0.1 max(app.bgObj.bgCh2)+0.1];
            plot(app.bg3UIAxes,app.bgObj.bgCh3,'Marker','.','LineStyle','--','Color', app.Ch3Color,'MarkerSize',8);
            title(app.bg3UIAxes,sprintf('Channel 3 Background, CtrlV = %.3f, Gain = %.3f',app.bgObj.CtrlV3, app.bgObj.bgGain3));
            app.bg3UIAxes.YLim = [min(app.bgObj.bgCh3)-0.1 max(app.bgObj.bgCh3)+0.1];
        end


        function save_LG_vs_mExp_error_fig(app)
            channel = app.channelDropDown.Value;
            switch channel
                case 'Channel 1'
                    plotObj = app.Ch1DataObj;
                    %                             c = app.Ch1Color;
                case 'Channel 2'
                    plotObj = app.Ch2DataObj;
                    %                             c = app.Ch2Color;
                case 'Channel 3'
                    plotObj = app.Ch3DataObj;
                    %                             c = app.Ch3Color;
            end
            res_Lg = get(plotObj,'wf_aligned')-get(plotObj,'fit');
            se_Lg = sum(res_Lg.^2)';
            res_Exp = plotObj.expDeconObj.spec_aligned-get(plotObj.expDeconObj,'fit');
            se_Exp = sum(res_Exp.^2)';
            f2 = figure('Position',[100 100 960 540]);
            ax = gca;
            scatter(ax,se_Lg,se_Exp, [],'magenta','filled','o');
            hold(ax,'on');
            plot(ax,[0 max(se_Exp)],[0 max(se_Exp)],'g--','LineWidth',1.5)
            hold(ax,'off');
            xlabel(ax,'Laguerre Square Error')
            ylabel(ax,'Muti-exponential Square Error')
            %                     axis(ax,'tight')
            title(ax,'Square Error')
            lg = {'Laguerre vs Exponential','y=x'};
            legend(ax,lg,'Location','northwest')
            ylim(ax,[quantile(se_Exp,0.05) quantile(se_Exp,0.95)])
            xlim(ax,[quantile(se_Lg,0.05) quantile(se_Lg,0.95)])
            filename = fullfile(app.dataFolder, app.saveFolderName,'LG_vs_mEXp_SquareError.fig');
            saveas(f2,filename)
            close(f2)
        end

        function loadSoftwareVersion(app)
            ini = IniConfig();
            ini.ReadFile('SoftwareVersion.ini');
            sections = ini.GetSections();
            [keys, ~] = ini.GetKeys(sections{1});
            temp = ini.GetValues(sections{1}, keys);
            MAJOR = temp{1};
            MINOR = temp{2};
            PATCH = temp{3};
            app.softwareVersion = [num2str(MAJOR) '.' num2str(MINOR) '.' num2str(PATCH)];
            app.UIFigure.Name = ['FLImBRUSH Data Processing Tool Version ' app.softwareVersion];
        end
    end





    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app)
            %             clear variables
            clc
            %             availableThreads = maxNumCompThreads;
            app.poolobj = gcp;
            if ~app.poolobj.Connected
                app.poolobj = parpool('threads');
            end
            loadSoftwareVersion(app)
            app.TabGroup.SelectedTab = app.fileSelectionTab;
            setAPDDetectorPath(app)
            loadDeConSettingIni(app)
            app.RemoveSaturationButton.Enable = 0;
            app.AverageDataButton.Enable = 0;
            app.RemoveBGButton.Enable = 0;
            app.TruncateButton.Enable = 0;
            app.DeconButton.Enable = 0;
            app.SaveDataButton.Enable = 0;
        end

        % Menu selected function: LoadRawDataMenu
        function LoadRawDataMenuSelected(app, event)
            app.irfAlignFlag = 0;
            [app.dataNameList, folderlist] = creatInputFileList(app.rootFolder);
            app.dataNameList = app.dataNameList';
            app.dataFolder = folderlist{1};
            app.dataFielesListBox.Items = app.dataNameList;
            app.dataFielesListBox.Value = app.dataNameList;
        end

        % Value changed function: bgNotNeededCheckBox
        function bgNotNeededCheckBoxValueChanged(app, event)
            value = app.bgNotNeededCheckBox.Value;
            if value
                app.bgLoadedLamp.Color = 'g';
                app.bgFileEditField.Value = 'Warning: No BG subtraction will be performed!';
                app.bgObj = []; % clear bg from memery
                cla(app.bg3UIAxes)
                cla(app.bg2UIAxes)
                cla(app.bg1UIAxes)
                app.bgLoadedFlag = 1;
            else
                app.bgLoadedLamp.Color = 'r';
                %                 app.bgFileEditField.Value = 'No Selection.';
            end
        end

        % Button pushed function: apd1LoadButton
        function apd1LoadButtonPushed(app, event)
            [filename, foldername] = uigetfile({'*.tdms;*.mat','APD Data Files (*.tdms,*.mat)'}, 'Please Select Channel 1 APD File', app.detectorDataFolder);
            figure(app.UIFigure)
            if filename
                app.detectorDataFolder = foldername;
                app.dataInfoObj.apd1Name = filename;
                app.dataInfoObj.apd1Folder = foldername;
                app.apd1Obj = apdClass(fullfile(foldername,filename));
                app.apd1LoadedLamp.Color = 'g';
                app.apd1LoadedFlag = 1;
                app.apd1FileEditField.Value = app.apd1Obj.filePath;
                creatFromPath(app.apd1Obj);
                app.apd1ViewButton.Enable = 'On'; % enable view button
            else
                app.apd1LoadedLamp.Color = 'r'; % reset indicator
                app.apd1FileEditField.Value = 'Please Select Channel 1 APD File'; % reset prompt
                app.apd1Obj = []; % clear irf field
                app.apd1LoadedFlag = 0; % clear flag
                app.apd1ViewButton.Enable = 'Off'; % enable view button
            end
        end

        % Button pushed function: loadDataButton
        function loadFileButtonPushed(app, event)
            [app.dataNameList, folderlist] = creatInputFileList(app.rootFolder);
            figure(app.UIFigure)
            if ~isempty(app.dataNameList)
                app.dataNameList = app.dataNameList';
                app.dataFolder = folderlist{1};
                app.dataFielesListBox.Items = app.dataNameList;
                app.dataFielesListBox.ItemsData = 1:length(app.dataNameList);
                app.dataFielesListBox.Value = app.fileSelection;
                app.DataFileSelectedDisp.Value = app.dataFielesListBox.Items{app.fileSelection};
                app.dataLoadedFlag = 1;
                app.clearButton.Enable = 'On'; % Enable clear button
                uiUpdate(app, 'init')
                clearAllPreviewAxis(app)
                clearAllFittingResultPlot(app) % clear fitting result figures
                app.processedDataloadedFlag = 0;
            else
                app.dataFielesListBox.Items = {'No Data File Selected!'};
                app.dataFielesListBox.Value = {'No Data File Selected!'};
                app.clearButton.Enable = 'Off'; % Disable clear button
            end
        end

        % Button pushed function: clearButton
        function clearButtonPushed(app, event)
            app.dataFielesListBox.Items = {'No Data file selected!'};
        end

        % Button pushed function: bgLoadButton
        function bgLoadButtonPushed(app, event)
            [filename_BG, foldername_BG] = uigetfile('*.tdms*', 'Please Select the Background tdms file', app.dataFolder);
            figure(app.UIFigure)
            %check whether file selection is valid
            if filename_BG % if valid, load bg
                app.dataInfoObj.bgFileName = filename_BG;
                app.dataInfoObj.bgFileFolder = foldername_BG;
                app.bgObj = backGround(fullfile(foldername_BG,filename_BG));
                loadBG(app.bgObj,app.upSampleFactor);
                app.bgObj.bgGain1 = interp1(app.apd1Obj.gainV,app.apd1Obj.apdGain,app.bgObj.CtrlV1,'spline');
                app.bgObj.bgGain2 = interp1(app.apd2Obj.gainV,app.apd2Obj.apdGain,app.bgObj.CtrlV2,'spline');
                app.bgObj.bgGain3 = interp1(app.apd3Obj.gainV,app.apd3Obj.apdGain,app.bgObj.CtrlV3,'spline');
                app.bgFileEditField.Value = app.bgObj.pathToBg;
                app.bgLoadedLamp.Color = 'g';
                app.bgLoadedFlag = 1;
                plotBG(app) % plot app
            else % if invalid, reset path indicator and lamp
                app.bgLoadedLamp.Color = 'r';
                app.bgFileEditField.Value = 'Please Select the Background tdms file!';
                app.bgObj = [];
                app.bgLoadedFlag = 0;
            end
        end

        % Selection change function: TabGroup
        function TabGroupSelectionChanged(app, event)
            %             selectedTab = app.TabGroup.SelectedTab;
            % check whether files are all configed.
            if (~app.apd1LoadedFlag || ~app.apd2LoadedFlag || ~app.apd3LoadedFlag || ~app.bgLoadedFlag || ~app.dataLoadedFlag) && ~app.processedDataloadedFlag
                uialert(app.UIFigure,'Files are not set, please make sure all necessary files are selected!','Invalid File');
                app.TabGroup.SelectedTab = app.fileSelectionTab;
            end

            if app.TabGroup.SelectedTab == app.DeConSettingTab
                uiUpdate(app,'init')
                clearBGAxis(app)
            end
            if (app.TabGroup.SelectedTab == app.fileSelectionTab) && ~app.processedDataloadedFlag
                plotBG(app)
            end
        end

        % Button pushed function: LoadDataButton
        function LoadDataButtonPushed(app, event)

            clearAllPreviewAxis(app) % clear all existing data plot
            clearAllFittingResultPlot(app) % clear all fitting result axis
            set(app.UIFigure,'pointer','watch') % set cursor to spinning
            drawnow
            app.dataInfoObj.dataFileNames = app.dataFielesListBox.Value;
            app.dataInfoObj.dataFolderName = app.dataFolder;
            app.FLImDataObj = FLImDataClass(fullfile(app.dataFolder,app.dataFielesListBox.Items{app.fileSelection}));
            loadFromfile(app.FLImDataObj)
            %--------------set laser rep rate field----------------------
            app.laserRepRateEditField.Value = app.FLImDataObj.laserRepRate;
            %--------------populate average selection----------------------
            avg = app.FLImDataObj.dataAvg;
            avg = double(avg);
            AvgSelection = cell(log2(avg)+1,1);
            for i = 1:log2(avg)+1
                AvgSelection{i} = num2str(avg/power(2,i-1));
            end
            app.DataAverageDropDown.Items = AvgSelection;
            app.DataAverageDropDown.Value = AvgSelection{1};
            %---------------------------------creat channel data object---------------------------------------------------------
            app.Ch1DataObj = ChannelDataAPD(app.FLImDataObj.ch1RawWF, app.FLImDataObj.V1, app.FLImDataObj.dataAvg, app.apd1Obj, app.TimeResolutionnsEditField.Value, app.bgObj.bgCh1, app.FLImDataObj.laserRepRate, app.upSampleFactor);
            app.Ch2DataObj = ChannelDataAPD(app.FLImDataObj.ch2RawWF, app.FLImDataObj.V2, app.FLImDataObj.dataAvg, app.apd2Obj, app.TimeResolutionnsEditField.Value, app.bgObj.bgCh2, app.FLImDataObj.laserRepRate, app.upSampleFactor);
            app.Ch3DataObj = ChannelDataAPD(app.FLImDataObj.ch3RawWF, app.FLImDataObj.V3, app.FLImDataObj.dataAvg, app.apd3Obj, app.TimeResolutionnsEditField.Value, app.bgObj.bgCh3, app.FLImDataObj.laserRepRate, app.upSampleFactor);
            removeDCData(app.Ch1DataObj);
            removeDCData(app.Ch2DataObj);
            removeDCData(app.Ch3DataObj);
            %----------------------plot raw data and marker-----------------
            plotData(app,'raw')
            %             plotDataIgnoreMarker(app)
            %             plotBGMarker(app)

            set(app.UIFigure,'pointer','arrow')
            drawnow

            uiUpdate(app, 'dataLoaded') % enable UI component
            app.SaveDataButton.Enable = 0;

        end

        % Value changed function: dataIgnoreLowEditField
        function dataIgnoreLowEditFieldValueChanged(app, event)
            plotDataIgnoreMarker(app)
        end

        % Value changed function: dataIgnoreHighEditField
        function dataIgnoreHighEditFieldValueChanged(app, event)
            plotDataIgnoreMarker(app)
        end

        % Value changed function: BgLowEditField
        function BgLowEditFieldValueChanged(app, event)
            plotBGMarker(app)
        end

        % Value changed function: BgHighEditField
        function BgHighEditFieldValueChanged(app, event)
            plotBGMarker(app)
        end

        % Button pushed function: RemoveBGButton
        function RemoveBGButtonPushed(app, event)

            set(app.UIFigure,'pointer','watch') % set cursor to spinning
            drawnow

            bgLow = app.BgLowEditField.Value;
            bgHigh = app.BgHighEditField.Value;
            removeDCBG(app.Ch1DataObj,bgLow,bgHigh,0.03);
            removeDCBG(app.Ch2DataObj,bgLow,bgHigh);
            removeDCBG(app.Ch3DataObj,bgLow,bgHigh);
            plotPreprocessedData(app)

            set(app.UIFigure,'pointer','arrow') % set cursor to arrow
            drawnow
            uiUpdate(app, 'BGRemoved')
        end

        % Callback function
        function TruncateDataButtonPushed(app, event)
            app.rawDataObj.peaks = [app.CH1EditField.Value app.CH2EditField.Value app.CH3EditField.Value app.CH4EditField.Value];
            truncatedData = truncateData(app.rawDataObj,app.channelWidthEditField.Value);
            app.ChDataArray{1} = ChannelData(truncatedData{1},app.irfObj.ch1Norm,app.TimeResolutionnsEditField.Value,app.DetectorBandwidthGHzEditField.Value,app.rawDataObj.noise,app.rawDataObj.gain);
            app.ChDataArray{2} = ChannelData(truncatedData{2},app.irfObj.ch2Norm,app.TimeResolutionnsEditField.Value,app.DetectorBandwidthGHzEditField.Value,app.rawDataObj.noise,app.rawDataObj.gain);
            app.ChDataArray{3} = ChannelData(truncatedData{3},app.irfObj.ch3Norm,app.TimeResolutionnsEditField.Value,app.DetectorBandwidthGHzEditField.Value,app.rawDataObj.noise,app.rawDataObj.gain);
            app.ChDataArray{4} = ChannelData(truncatedData{4},app.irfObj.ch4Norm,app.TimeResolutionnsEditField.Value,app.DetectorBandwidthGHzEditField.Value,app.rawDataObj.noise,app.rawDataObj.gain);
            app.TabGroup.SelectedTab = app.channelDataTab;
            plotChannelData(app, app.UIAxes_Ch1, app.ChDataArray{1});
            plotChannelData(app, app.UIAxes_Ch2, app.ChDataArray{2});
            plotChannelData(app, app.UIAxes_Ch3, app.ChDataArray{3});
            plotChannelData(app, app.UIAxes_Ch4, app.ChDataArray{4});
        end

        % Callback function
        function AlignIRFButtonPushed(app, event)
            % get channel selection
            app.deconChannelSelection = [app.channel1CheckBox.Value app.channel2CheckBox.Value*2 app.channel3CheckBox.Value*3 app.channel4CheckBox.Value*4];
            app.deconChannelSelection(app.deconChannelSelection==0)=[];
            app.ChLaguerreArray = cell(4,1); % clear cell array that stores the obj
            app.irfAlignFlag = 0; %clear irf alignment flag
            for i=1:4 % assign Laguerre obj to cell array
                app.ChLaguerreArray{i} = LaguerreModel(app.ChDataArray{i});
            end
            % align irf
            for i=app.deconChannelSelection
                iIRF_align(app.ChLaguerreArray{i});
            end
            % plot irf alignment result

            plotChannelData(app, app.UIAxes_Ch1, app.ChLaguerreArray{1}.channeldata);
            plotChannelData(app, app.UIAxes_Ch2, app.ChLaguerreArray{2}.channeldata);
            plotChannelData(app, app.UIAxes_Ch3, app.ChLaguerreArray{3}.channeldata);
            plotChannelData(app, app.UIAxes_Ch4, app.ChLaguerreArray{4}.channeldata);

            % set irf alignment flag to 1.
            app.irfAlignFlag = 1;
        end

        % Callback function
        function RunDeconButtonPushed(app, event)
            % set decon indicator and string
            app.deConStatusLamp.Color = 'r';
            app.deConStatusLabel.Text = 'Busy';
            drawnow % update UI
            % run decon for selected channel
            for i =app.deconChannelSelection
                if ~app.irfAlignFlag
                    iIRF_align(app.ChLaguerreArray{i});
                end
                app.irfAlignFlag = 1;
                dataThresholding(app.ChLaguerreArray{i},[app.LowThresholdEditField.Value app.HighThresholdEditField.Value]);
                estimate_laguerre(app.ChLaguerreArray{i});
            end
            % reset decon indicator and text
            app.deConStatusLamp.Color = 'g';
            app.deConStatusLabel.Text = 'Idle';
            % display decon done in command window
            disp('Deconvolution Done!')
            % switch to decon fitting overview tab
            app.TabGroup.SelectedTab = app.FittingResultTab;
            % set decon fitting channel selector items
            Items_temp = cell(size(app.deconChannelSelection));
            for i = 1:length(app.deconChannelSelection)
                Items_temp{i} = num2str(app.deconChannelSelection(i));
            end
            app.channelDropDown.Items = Items_temp;
            % plot deconvolusion result on tab
            plotDeconResult(app)
            % creat decon result object that stores all data and reconstruct all images.
            app.deconResutlObj = deconResult(app.rawDataObj,app.ChLaguerreArray);
            %             plotImg(app)
        end

        % Value changed function: LineSlider
        function LineSliderValueChanged(app, event)
            % to validate integer value
            value = app.LineSlider.Value;
            value  = round(value);
            app.LineSlider.Value = value;
            app.PointEditField.Value = value;
            plotDeconResult(app)
            option = app.PlotDropDown.Value;
            if strcmp(option,'Raw Waveforms')
                plotTrace(app, option)
            end
        end

        % Value changed function: channelDropDown
        function channelDropDownValueChanged(app, event)
            plotDeconResult(app)
            option = app.PlotDropDown.Value;
            plotTrace(app,option)
            %             plotImg(app)
        end

        % Value changed function: xMaxEditField, xMinEditField
        function xMinEditFieldValueChanged(app, event)
            plotDeconResult(app)

        end

        % Callback function
        function CH1EditFieldValueChanged(app, event)
            % plop new peak location
            plotPeakLocations(app, [app.CH1EditField.Value app.CH2EditField.Value app.CH3EditField.Value app.CH4EditField.Value]);
        end

        % Button pushed function: apd2LoadButton
        function apd2LoadButtonPushed(app, event)
            [filename, foldername] = uigetfile({'*.tdms;*.mat','APD Data Files (*.tdms,*.mat)'}, 'Please Select Channel 2 APD File', app.detectorDataFolder);
            figure(app.UIFigure)
            if filename
                app.detectorDataFolder = foldername;
                app.dataInfoObj.apd2Name = filename;
                app.dataInfoObj.apd2Folder = foldername;
                app.apd2Obj = apdClass(fullfile(foldername,filename));
                app.apd2LoadedLamp.Color = 'g';
                app.apd2LoadedFlag = 1;
                app.apd2FileEditField.Value = app.apd2Obj.filePath;
                creatFromPath(app.apd2Obj);
                app.apd2ViewButton.Enable = 'On'; % enable view button
            else
                app.apd2LoadedLamp.Color = 'r'; % reset indicator
                app.apd2FileEditField.Value = 'Please Select Channel 2 APD File'; % reset prompt
                app.apd2Obj = []; % clear irf field
                app.apd2LoadedFlag = 0; % clear flag
                app.apd2ViewButton.Enable = 'Off'; % Disable view button
            end
        end

        % Button pushed function: apd3LoadButton
        function apd3LoadButtonPushed(app, event)
            [filename, foldername] = uigetfile({'*.tdms;*.mat','APD Data Files (*.tdms;*.mat)'}, 'Please Select Channel 3 APD File', app.detectorDataFolder);
            figure(app.UIFigure)
            if filename
                app.detectorDataFolder = foldername;
                app.dataInfoObj.apd3Name = filename;
                app.dataInfoObj.apd3Folder = foldername;
                app.apd3Obj = apdClass(fullfile(foldername,filename));
                app.apd3LoadedLamp.Color = 'g';
                app.apd3LoadedFlag = 1;
                app.apd3FileEditField.Value = app.apd3Obj.filePath;
                creatFromPath(app.apd3Obj);
                app.apd3ViewButton.Enable = 'On'; % enable view button
            else
                app.apd3LoadedLamp.Color = 'r'; % reset indicator
                app.apd3FileEditField.Value = 'Please Select Channel 3 APD File'; % reset prompt
                app.apd3Obj = []; % clear irf field
                app.apd3LoadedFlag = 0; % clear flag
                app.apd3ViewButton.Enable = 'Off'; % Disable view button
            end
        end

        % Button pushed function: caliLoadButton
        function caliLoadButtonPushed(app, event)
            [filename, foldername] = uigetfile('*.tdms', 'Please Select System Calibration File', app.dataFolder);
            figure(app.UIFigure)
            if filename
                app.dataInfoObj.caliFileName = filename;
                app.dataInfoObj.caliFileFolder = foldername;
                app.caliObj = sysCaliDataClass(fullfile(foldername,filename));
                app.caliFileEditField.Value = app.caliObj.filePath;
                app.caliLoadedLamp.Color = 'g';
                load(app.caliObj);
            else
                app.caliLoadedLamp.Color = 'r'; % reset indicator
                app.caliFileEditField.Value = 'Please Select System Calibration File'; % reset prompt
                app.caliObj = []; % clear irf field
                app.caliLoadedFlag = 0; % clear flag
            end

        end

        % Button pushed function: apd1ViewButton
        function apd1ViewButtonPushed(app, event)
            app.apd1ViewButton.Enable = 'Off';
            app.apd1DialogApp = appApdViewInfo(app,app.apd1Obj,1);

        end

        % Close request function: UIFigure
        function MainAppCloseRequest(app, event)
            saveAPDDetectorPath(app);
            saveDeconSettingIni(app);
            delete(app.apd1DialogApp);
            delete(app.apd2DialogApp);
            delete(app.apd3DialogApp);
            %             delete(app.poolobj); % delete pool object
            delete(app);

        end

        % Button pushed function: apd2ViewButton
        function apd2ViewButtonPushed(app, event)
            app.apd2ViewButton.Enable = 'Off';
            app.apd2DialogApp = appApdViewInfo(app,app.apd2Obj,2);
        end

        % Button pushed function: apd3ViewButton
        function apd3ViewButtonPushed(app, event)
            app.apd3ViewButton.Enable = 'Off';
            app.apd3DialogApp = appApdViewInfo(app,app.apd3Obj,3);
        end

        % Button pushed function: RemoveSaturationButton
        function RemoveSaturationButtonPushed(app, event)
            set(app.UIFigure,'pointer','watch') % set cursor to spinning
            drawnow

            [totalNumOfWf, goodNumOfWf, badNumOfWf] = removeSaturation(app.Ch1DataObj,app.LowThresholdEditField.Value,app.HighThresholdEditField.Value);
            msg1 = sprintf('Total number of waveforms: %d.\nWaveforms left: %d.\nWaveforms removed: %d.',...
                totalNumOfWf,goodNumOfWf,badNumOfWf);
            [totalNumOfWf, goodNumOfWf, badNumOfWf] = removeSaturation(app.Ch2DataObj,app.LowThresholdEditField.Value,app.HighThresholdEditField.Value);
            msg2 = sprintf('Total number of waveforms: %d.\nWaveforms left: %d.\nWaveforms removed: %d.',...
                totalNumOfWf,goodNumOfWf,badNumOfWf);
            [totalNumOfWf, goodNumOfWf, badNumOfWf] = removeSaturation(app.Ch3DataObj,app.LowThresholdEditField.Value,app.HighThresholdEditField.Value);
            msg3 = sprintf('Total number of waveforms: %d.\nWaveforms left: %d.\nWaveforms removed: %d.',...
                totalNumOfWf,goodNumOfWf,badNumOfWf);
            msgbox({'Channel 1:';msg1;'';'Channel 2:';msg2;'';'Channel 3:';msg3},'Data info','help')
            plotData(app,'filtered')
            %             plotDataIgnoreMarker(app)
            %             plotBGMarker(app)

            app.AverageDataButton.Enable = 1;

            set(app.UIFigure,'pointer','arrow') % set cursor to arrow
            drawnow
            uiUpdate(app, 'removedSat')

        end

        % Button pushed function: AverageDataButton
        function AverageDataButtonPushed(app, event)

            set(app.UIFigure,'pointer','watch') % set cursor to spinning
            drawnow

            avg = app.DataAverageDropDown.Value;
            avg = str2double(avg);
            avgData(app.Ch1DataObj, avg)
            avgData(app.Ch2DataObj, avg)
            avgData(app.Ch3DataObj, avg)
            plotData(app,'averaged')
            %             plotDataIgnoreMarker(app)
            plotBGMarker(app)

            app.RemoveBGButton.Enable = 1;

            set(app.UIFigure,'pointer','arrow') % set cursor to arrow
            drawnow
            uiUpdate(app, 'dataAveraged')

        end

        % Button pushed function: TruncateButton
        function TruncateButtonPushed(app, event)
            %             laguerreOrder = app.LaguerreOrderEditField.Value;
            %             laguerreAlpha = app.Ch1AlphaValueDropDown.Value;
            %             minTrunLength = findTruncationL(laguerreOrder,laguerreAlpha);
            %             minTimeWindow = app.Ch1DataObj.dtUp*minTrunLength;
            %             if app.channelWidthEditField.Value<minTimeWindow
            %                 message1 = sprintf('Truncation time windows too short for alpha value of %.4f. \n',laguerreAlpha);
            %                 message2 = sprintf('Minimum truncation windwos of %.2f will be automatic set.\n',minTimeWindow);
            %                 message3 = 'if data is not long enough, 0 will be added in the tail.';
            % %                 app.channelWidthEditField.Value = minTimeWindow; %do not overwrite window
            %                 f = msgbox([message1 message2 message3]);
            %             end

            set(app.UIFigure,'pointer','watch') % set cursor to spinning
            drawnow

            %             upSampleData(app.Ch1DataObj);
            truncateData(app.Ch1DataObj,app.channelWidthEditField.Value);
            %             upSampleData(app.Ch2DataObj);
            truncateData(app.Ch2DataObj,app.channelWidthEditField.Value);
            %             upSampleData(app.Ch3DataObj);
            truncateData(app.Ch3DataObj,app.channelWidthEditField.Value);
            plotTruncatedData(app)

            plotDataIgnoreMarker(app)

            set(app.UIFigure,'pointer','arrow') % set cursor to arrow
            drawnow
            uiUpdate(app, 'truncated')
        end

        % Button pushed function: DeconButton
        function DeconButtonPushed(app, event)

            laguerreOrder = app.LaguerreOrderEditField.Value;
            alpha1 = app.Ch1AlphaValueDropDown.Value;
            alpha2 = app.Ch2AlphaValueDropDown.Value;
            alpha3 = app.Ch3AlphaValueDropDown.Value;
            alpha1 = str2double(alpha1);
            alpha2 = str2double(alpha2);
            alpha3 = str2double(alpha3);

            %             alphaMax = alpha_up(app.channelWidthEditField.Value/0.08,laguerreOrder)+0.1;
            alphaMax = 1;
            if any(alphaMax<[alpha1 alpha2 alpha3])
                msg1 = sprintf('Maxmum allowed alpha value for your truncation length is %.3f \n', alphaMax);
                msg2 = 'Please reselect your alpha value!\n';
                msgbox([msg1 msg2],'modal');
            else

                exclude = app.dataIgnoreLowEditField.Value:app.dataIgnoreHighEditField.Value;
                if any(exclude < 1)
                    exclude = [];
                end

                set(app.UIFigure,'pointer','watch') % set cursor to spinning
                drawnow

                f = waitbar(0,'Running Laguerre channel 1...');
                runDeconLG(app.Ch1DataObj,exclude,laguerreOrder,alpha1);

                waitbar(0.33,f,'Running Laguerre channel 2...');
                runDeconLG(app.Ch2DataObj,exclude,laguerreOrder,alpha2);

                waitbar(0.66,f,'Running Laguerre channel 3...');
                runDeconLG(app.Ch3DataObj,exclude,laguerreOrder,alpha3);
                waitbar(1,f,'Laguerre fit finished...');
                
                app.expDeconFlag = app.RunEXPCheckBox.Value;

                if app.expDeconFlag
                    expOrder = app.exponentialsDropDown.Value;
                    expOrder = str2double(expOrder);
                    f1 = waitbar(0,'Running exponential channel 1...');
                    runDeconExp(app.Ch1DataObj,expOrder,[],[],[],exclude);
                    waitbar(0.33,f1,'Running exponential channel 2...');
                    runDeconExp(app.Ch2DataObj,expOrder,[],[],[],exclude);
                    waitbar(0.66,f1,'Running exponential channel 3...');
                    runDeconExp(app.Ch3DataObj,expOrder,[],[],[],exclude);
                    waitbar(1,f1,'exponential fit finished...');
                end

                %------------------run phasor for harmonic of 1 to 4-----------------------
                for H = 1:4 % run phasor for harmonic of 1 to 4
                    runPhasor(app.Ch1DataObj, H);
                    runPhasor(app.Ch2DataObj, H);
                    runPhasor(app.Ch3DataObj, H);
                end

                CH1WFGainCorr = get(app.Ch1DataObj,'GainCorrectedWF');
                CH2WFGainCorr = get(app.Ch2DataObj,'GainCorrectedWF');
                CH3WFGainCorr = get(app.Ch3DataObj,'GainCorrectedWF');
                %-------------------------------------------compute extended phasor-----------------------------------------------------
                WFExtended = [CH1WFGainCorr;CH2WFGainCorr;CH3WFGainCorr];
                numOfWF = size(WFExtended,2);
                EOP = zeros(numOfWF,1);
                for i = 1:numOfWF
                    EOP(i) = ComputePhasor(WFExtended(:,i),0,1,0);
                end
                app.EOP_H1G = real(EOP);
                app.EOP_H1S = imag(EOP);
                %-------------------------------------------compute spectral phasor--------------------------------------------
                Spectrum = zeros(numOfWF,3);
                Spectrum(:,1) = sum(CH1WFGainCorr)';
                Spectrum(:,2) = sum(CH2WFGainCorr)';
                Spectrum(:,3) = sum(CH3WFGainCorr)';
                SP = zeros(numOfWF,1);
                for i = 1: numOfWF
                    SP(i) = (ComputePhasor(Spectrum(i,:),0,1,0)+(1+1i))/2;
                end
                app.SP_G = real(SP);
                app.SP_S = imag(SP);

                waitbar(1,f,'Finished');

                close(f)
                delete(f)
                % close(f1)
                % delete(f1)

                set(app.UIFigure,'pointer','arrow') % set cursor to arrow
                drawnow


                app.TabGroup.SelectedTab = app.FittingResultTab;
                clearAllPreviewAxis(app)
                app.LineSlider.Value = 1;
                updateChannelDropDown(app)
                plotDeconResult(app);
                option = app.PlotDropDown.Value;
                plotTrace(app,option)

                app.SaveDataButton.Enable = 1;
            end



        end

        % Button pushed function: SaveDataButton
        function SaveDataButtonPushed(app, event)

            set(app.UIFigure,'pointer','watch') % set cursor to spining
            drawnow

            c = clock; % get current time and date
            app.saveFolderName = ['Decon ' num2str(c(1)) '-' num2str(c(2)) '-' num2str(c(3)) ' ' num2str(c(4)) '-' num2str(c(5))];
            mkdir(app.dataFolder,app.saveFolderName);
            [~,saveFileName,~] = fileparts(app.DataFileSelectedDisp.Value);
            saveFileName = [saveFileName '.mat'];
            saveFileFullPath = fullfile(app.dataFolder,app.saveFolderName,saveFileName);

            dataInfoObj = app.dataInfoObj;

            Ch1DataObj = app.Ch1DataObj;
            Ch2DataObj = app.Ch2DataObj;
            Ch3DataObj = app.Ch3DataObj;
            EOP_H1G = app.EOP_H1G;
            EOP_H1S = app.EOP_H1S;
            SP_G = app.SP_G;
            SP_S = app.SP_S;
            calibrationObj =app.caliObj;

            save(saveFileFullPath, 'dataInfoObj','Ch1DataObj','Ch2DataObj','Ch3DataObj','calibrationObj','EOP_H1G','EOP_H1S','SP_G','SP_S','-v7.3');
            saveDeconLite(Ch1DataObj,Ch2DataObj,Ch3DataObj,saveFileFullPath)
            %save Laguerre vs multi-exp square error
            %             save_LG_vs_mExp_lifetime_fig(app)
            %             save_LG_vs_mExp_error_fig(app)

            set(app.UIFigure,'pointer','arrow') % set cursor to arrow
            drawnow

            msgbox('Data Successfully Saved!')
        end

        % Value changed function: PlotDropDown
        function PlotDropDownValueChanged(app, event)
            option = app.PlotDropDown.Value;
            plotTrace(app,option)
        end

        % Button pushed function: UpsampleButton
        function UpsampleButtonPushed(app, event)
            set(app.UIFigure,'pointer','watch') % set cursor to spinning
            drawnow
            upSampleData(app.Ch1DataObj);
            upSampleData(app.Ch2DataObj);
            upSampleData(app.Ch3DataObj);
            alignWF_CFD(app.Ch1DataObj, 1, (80:180)*app.upSampleFactor)
            alignWF_CFD(app.Ch2DataObj, 0.5, (180:size(app.Ch2DataObj.rawData,1))*app.upSampleFactor)
            alignWF_CFD(app.Ch3DataObj, 0.5, (180:size(app.Ch2DataObj.rawData,1))*app.upSampleFactor)
            plotUpsampledData(app)
            uiUpdate(app,'Upsampled')
            set(app.UIFigure,'pointer','arrow') % set cursor to spinning
            drawnow
        end

        % Value changed function: dataFielesListBox
        function dataFielesListBoxValueChanged(app, event)
            app.fileSelection = app.dataFielesListBox.Value;
            app.DataFileSelectedDisp.Value = app.dataFielesListBox.Items{app.fileSelection};
        end

        % Value changed function: PointEditField
        function PointEditFieldValueChanged(app, event)
            value = app.PointEditField.Value;
            value = round(value);
            app.LineSlider.Value = value;
            plotDeconResult(app)
            option = app.PlotDropDown.Value;
            if strcmp(option,'Raw Waveforms')
                plotTrace(app, option)
            end

        end

        % Value changed function: channelWidthEditField
        function channelWidthEditFieldValueChanged(app, event)
            plotPreprocessedData(app)
        end

        % Menu selected function: LoadProcessedDataMenu
        function LoadProcessedDataMenuSelected(app, event)
            [file,path] = uigetfile(app.rootFolder,'*.mat','MultiSelect',"off",'Please select deconvolved .mat file.');
            temp = load(fullfile(path,file));
            app.Ch1DataObj = temp.Ch1DataObj;
            app.Ch2DataObj = temp.Ch1DataObj;
            app.Ch3DataObj = temp.Ch2DataObj;
            app.dataInfoObj = temp.dataInfoObj;
            app.caliObj = temp.calibrationObj;
            app.processedDataloadedFlag = 1;
            app.TabGroup.SelectedTab = app.FittingResultTab;
            clearAllPreviewAxis(app)
            app.LineSlider.Value = 1;
            updateChannelDropDown(app)
            plotDeconResult(app);
            option = app.PlotDropDown.Value;
            plotTrace(app,option)
        end

        % Menu selected function: SetRootPathMenu
        function SetRootPathMenuSelected(app, event)
            app.rootFolder = uigetdir;
        end

        % Value changed function: RunEXPCheckBox
        function RunEXPCheckBoxValueChanged(app, event)
            value = app.RunEXPCheckBox.Value;
            if value
                app.exponentialsDropDown.Enable = 1;
            else
                app.exponentialsDropDown.Enable = 0;
            end
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [100 100 1247 727];
            app.UIFigure.Name = 'Error';
            app.UIFigure.Icon = 'Icon.png';
            app.UIFigure.CloseRequestFcn = createCallbackFcn(app, @MainAppCloseRequest, true);

            % Create MainApp
            app.MainApp = uimenu(app.UIFigure);
            app.MainApp.Text = 'Data';

            % Create LoadRawDataMenu
            app.LoadRawDataMenu = uimenu(app.MainApp);
            app.LoadRawDataMenu.MenuSelectedFcn = createCallbackFcn(app, @LoadRawDataMenuSelected, true);
            app.LoadRawDataMenu.Text = 'Load Raw Data';

            % Create LoadProcessedDataMenu
            app.LoadProcessedDataMenu = uimenu(app.MainApp);
            app.LoadProcessedDataMenu.MenuSelectedFcn = createCallbackFcn(app, @LoadProcessedDataMenuSelected, true);
            app.LoadProcessedDataMenu.Text = 'Load Processed Data';

            % Create SetRootPathMenu
            app.SetRootPathMenu = uimenu(app.MainApp);
            app.SetRootPathMenu.MenuSelectedFcn = createCallbackFcn(app, @SetRootPathMenuSelected, true);
            app.SetRootPathMenu.Text = 'Set Root Path';

            % Create ExportMenu
            app.ExportMenu = uimenu(app.UIFigure);
            app.ExportMenu.Text = 'Export';

            % Create TabGroup
            app.TabGroup = uitabgroup(app.UIFigure);
            app.TabGroup.SelectionChangedFcn = createCallbackFcn(app, @TabGroupSelectionChanged, true);
            app.TabGroup.Position = [6 5 1240 720];

            % Create fileSelectionTab
            app.fileSelectionTab = uitab(app.TabGroup);
            app.fileSelectionTab.Title = 'File Selection';

            % Create GridLayout
            app.GridLayout = uigridlayout(app.fileSelectionTab);
            app.GridLayout.ColumnWidth = {60, 131, 128, 63, '29.09x', '0x', 40, 83, 23, '3x'};
            app.GridLayout.RowHeight = {22, 22, 22, 22, 22, '1x', '1x', '1x'};
            app.GridLayout.RowSpacing = 6.55555555555556;
            app.GridLayout.Padding = [10 6.55555555555556 10 6.55555555555556];

            % Create bg1UIAxes
            app.bg1UIAxes = uiaxes(app.GridLayout);
            title(app.bg1UIAxes, 'Ch1 background')
            xlabel(app.bg1UIAxes, 'Points')
            ylabel(app.bg1UIAxes, 'Voltage(V)')
            app.bg1UIAxes.XTickLabelRotation = 0;
            app.bg1UIAxes.YTickLabelRotation = 0;
            app.bg1UIAxes.ZTickLabelRotation = 0;
            app.bg1UIAxes.Layout.Row = 6;
            app.bg1UIAxes.Layout.Column = [4 10];

            % Create bg2UIAxes
            app.bg2UIAxes = uiaxes(app.GridLayout);
            title(app.bg2UIAxes, 'Ch2 background')
            xlabel(app.bg2UIAxes, 'Points')
            ylabel(app.bg2UIAxes, 'Voltage(V)')
            app.bg2UIAxes.XTickLabelRotation = 0;
            app.bg2UIAxes.YTickLabelRotation = 0;
            app.bg2UIAxes.ZTickLabelRotation = 0;
            app.bg2UIAxes.Layout.Row = 7;
            app.bg2UIAxes.Layout.Column = [4 10];

            % Create bg3UIAxes
            app.bg3UIAxes = uiaxes(app.GridLayout);
            title(app.bg3UIAxes, 'Ch3 background')
            xlabel(app.bg3UIAxes, 'Points')
            ylabel(app.bg3UIAxes, 'Voltage(V)')
            app.bg3UIAxes.XTickLabelRotation = 0;
            app.bg3UIAxes.YTickLabelRotation = 0;
            app.bg3UIAxes.ZTickLabelRotation = 0;
            app.bg3UIAxes.Layout.Row = 8;
            app.bg3UIAxes.Layout.Column = [4 10];

            % Create caliViewButton
            app.caliViewButton = uibutton(app.GridLayout, 'push');
            app.caliViewButton.BusyAction = 'cancel';
            app.caliViewButton.Enable = 'off';
            app.caliViewButton.Layout.Row = 4;
            app.caliViewButton.Layout.Column = [9 10];
            app.caliViewButton.Text = 'View';

            % Create caliLoadButton
            app.caliLoadButton = uibutton(app.GridLayout, 'push');
            app.caliLoadButton.ButtonPushedFcn = createCallbackFcn(app, @caliLoadButtonPushed, true);
            app.caliLoadButton.BusyAction = 'cancel';
            app.caliLoadButton.Layout.Row = 4;
            app.caliLoadButton.Layout.Column = 8;
            app.caliLoadButton.Text = 'Load';

            % Create caliLoadedLamp
            app.caliLoadedLamp = uilamp(app.GridLayout);
            app.caliLoadedLamp.Layout.Row = 4;
            app.caliLoadedLamp.Layout.Column = 7;
            app.caliLoadedLamp.Color = [1 0 0];

            % Create apd3ViewButton
            app.apd3ViewButton = uibutton(app.GridLayout, 'push');
            app.apd3ViewButton.ButtonPushedFcn = createCallbackFcn(app, @apd3ViewButtonPushed, true);
            app.apd3ViewButton.BusyAction = 'cancel';
            app.apd3ViewButton.Enable = 'off';
            app.apd3ViewButton.Layout.Row = 3;
            app.apd3ViewButton.Layout.Column = [9 10];
            app.apd3ViewButton.Text = 'View';

            % Create apd2ViewButton
            app.apd2ViewButton = uibutton(app.GridLayout, 'push');
            app.apd2ViewButton.ButtonPushedFcn = createCallbackFcn(app, @apd2ViewButtonPushed, true);
            app.apd2ViewButton.BusyAction = 'cancel';
            app.apd2ViewButton.Enable = 'off';
            app.apd2ViewButton.Layout.Row = 2;
            app.apd2ViewButton.Layout.Column = [9 10];
            app.apd2ViewButton.Text = 'View';

            % Create apd1ViewButton
            app.apd1ViewButton = uibutton(app.GridLayout, 'push');
            app.apd1ViewButton.ButtonPushedFcn = createCallbackFcn(app, @apd1ViewButtonPushed, true);
            app.apd1ViewButton.BusyAction = 'cancel';
            app.apd1ViewButton.Enable = 'off';
            app.apd1ViewButton.Layout.Row = 1;
            app.apd1ViewButton.Layout.Column = [9 10];
            app.apd1ViewButton.Text = 'View';

            % Create apd3LoadButton
            app.apd3LoadButton = uibutton(app.GridLayout, 'push');
            app.apd3LoadButton.ButtonPushedFcn = createCallbackFcn(app, @apd3LoadButtonPushed, true);
            app.apd3LoadButton.BusyAction = 'cancel';
            app.apd3LoadButton.Layout.Row = 3;
            app.apd3LoadButton.Layout.Column = 8;
            app.apd3LoadButton.Text = 'Load';

            % Create apd3LoadedLamp
            app.apd3LoadedLamp = uilamp(app.GridLayout);
            app.apd3LoadedLamp.Layout.Row = 3;
            app.apd3LoadedLamp.Layout.Column = 7;
            app.apd3LoadedLamp.Color = [1 0 0];

            % Create apd2LoadButton
            app.apd2LoadButton = uibutton(app.GridLayout, 'push');
            app.apd2LoadButton.ButtonPushedFcn = createCallbackFcn(app, @apd2LoadButtonPushed, true);
            app.apd2LoadButton.BusyAction = 'cancel';
            app.apd2LoadButton.Layout.Row = 2;
            app.apd2LoadButton.Layout.Column = 8;
            app.apd2LoadButton.Text = 'Load';

            % Create apd2LoadedLamp
            app.apd2LoadedLamp = uilamp(app.GridLayout);
            app.apd2LoadedLamp.Layout.Row = 2;
            app.apd2LoadedLamp.Layout.Column = 7;
            app.apd2LoadedLamp.Color = [1 0 0];

            % Create clearButton
            app.clearButton = uibutton(app.GridLayout, 'push');
            app.clearButton.ButtonPushedFcn = createCallbackFcn(app, @clearButtonPushed, true);
            app.clearButton.Enable = 'off';
            app.clearButton.Layout.Row = 1;
            app.clearButton.Layout.Column = 3;
            app.clearButton.Text = 'Clear';

            % Create loadDataButton
            app.loadDataButton = uibutton(app.GridLayout, 'push');
            app.loadDataButton.ButtonPushedFcn = createCallbackFcn(app, @loadFileButtonPushed, true);
            app.loadDataButton.BusyAction = 'cancel';
            app.loadDataButton.Layout.Row = 1;
            app.loadDataButton.Layout.Column = 2;
            app.loadDataButton.Text = 'Select';

            % Create bgNotNeededCheckBox
            app.bgNotNeededCheckBox = uicheckbox(app.GridLayout);
            app.bgNotNeededCheckBox.ValueChangedFcn = createCallbackFcn(app, @bgNotNeededCheckBoxValueChanged, true);
            app.bgNotNeededCheckBox.Text = 'BG not needed ';
            app.bgNotNeededCheckBox.Layout.Row = 5;
            app.bgNotNeededCheckBox.Layout.Column = [9 10];

            % Create bgLoadButton
            app.bgLoadButton = uibutton(app.GridLayout, 'push');
            app.bgLoadButton.ButtonPushedFcn = createCallbackFcn(app, @bgLoadButtonPushed, true);
            app.bgLoadButton.BusyAction = 'cancel';
            app.bgLoadButton.Layout.Row = 5;
            app.bgLoadButton.Layout.Column = 8;
            app.bgLoadButton.Text = 'Load';

            % Create apd1LoadButton
            app.apd1LoadButton = uibutton(app.GridLayout, 'push');
            app.apd1LoadButton.ButtonPushedFcn = createCallbackFcn(app, @apd1LoadButtonPushed, true);
            app.apd1LoadButton.BusyAction = 'cancel';
            app.apd1LoadButton.Layout.Row = 1;
            app.apd1LoadButton.Layout.Column = 8;
            app.apd1LoadButton.Text = 'Load';

            % Create bgLoadedLamp
            app.bgLoadedLamp = uilamp(app.GridLayout);
            app.bgLoadedLamp.Layout.Row = 5;
            app.bgLoadedLamp.Layout.Column = 7;
            app.bgLoadedLamp.Color = [1 0 0];

            % Create apd1FileEditField
            app.apd1FileEditField = uieditfield(app.GridLayout, 'text');
            app.apd1FileEditField.Editable = 'off';
            app.apd1FileEditField.HorizontalAlignment = 'right';
            app.apd1FileEditField.Layout.Row = 1;
            app.apd1FileEditField.Layout.Column = [5 6];
            app.apd1FileEditField.Value = 'Please Select Channel 1 APD File!';

            % Create APD1Label
            app.APD1Label = uilabel(app.GridLayout);
            app.APD1Label.HorizontalAlignment = 'right';
            app.APD1Label.Layout.Row = 1;
            app.APD1Label.Layout.Column = 4;
            app.APD1Label.Text = 'APD1';

            % Create bgFileEditField
            app.bgFileEditField = uieditfield(app.GridLayout, 'text');
            app.bgFileEditField.Editable = 'off';
            app.bgFileEditField.HorizontalAlignment = 'center';
            app.bgFileEditField.Layout.Row = 5;
            app.bgFileEditField.Layout.Column = [5 6];
            app.bgFileEditField.Value = 'Please Select Background File!';

            % Create BGfileEditFieldLabel
            app.BGfileEditFieldLabel = uilabel(app.GridLayout);
            app.BGfileEditFieldLabel.HorizontalAlignment = 'right';
            app.BGfileEditFieldLabel.Layout.Row = 5;
            app.BGfileEditFieldLabel.Layout.Column = 4;
            app.BGfileEditFieldLabel.Text = 'BG file:';

            % Create apd2FileEditField
            app.apd2FileEditField = uieditfield(app.GridLayout, 'text');
            app.apd2FileEditField.Editable = 'off';
            app.apd2FileEditField.HorizontalAlignment = 'right';
            app.apd2FileEditField.Layout.Row = 2;
            app.apd2FileEditField.Layout.Column = [5 6];
            app.apd2FileEditField.Value = 'Please Select Channel 2 APD File!';

            % Create APD2Label
            app.APD2Label = uilabel(app.GridLayout);
            app.APD2Label.HorizontalAlignment = 'right';
            app.APD2Label.Layout.Row = 2;
            app.APD2Label.Layout.Column = 4;
            app.APD2Label.Text = 'APD2';

            % Create apd3FileEditField
            app.apd3FileEditField = uieditfield(app.GridLayout, 'text');
            app.apd3FileEditField.Editable = 'off';
            app.apd3FileEditField.HorizontalAlignment = 'right';
            app.apd3FileEditField.Layout.Row = 3;
            app.apd3FileEditField.Layout.Column = [5 6];
            app.apd3FileEditField.Value = 'Please Select Channel 3 APD File!';

            % Create APD3Label
            app.APD3Label = uilabel(app.GridLayout);
            app.APD3Label.HorizontalAlignment = 'right';
            app.APD3Label.Layout.Row = 3;
            app.APD3Label.Layout.Column = 4;
            app.APD3Label.Text = 'APD3';

            % Create caliFileEditField
            app.caliFileEditField = uieditfield(app.GridLayout, 'text');
            app.caliFileEditField.Editable = 'off';
            app.caliFileEditField.HorizontalAlignment = 'center';
            app.caliFileEditField.Layout.Row = 4;
            app.caliFileEditField.Layout.Column = [5 6];
            app.caliFileEditField.Value = 'Please Select System Intensity Calibration!';

            % Create CalibrationLabel
            app.CalibrationLabel = uilabel(app.GridLayout);
            app.CalibrationLabel.HorizontalAlignment = 'right';
            app.CalibrationLabel.Layout.Row = 4;
            app.CalibrationLabel.Layout.Column = 4;
            app.CalibrationLabel.Text = 'Calibration';

            % Create apd1LoadedLamp
            app.apd1LoadedLamp = uilamp(app.GridLayout);
            app.apd1LoadedLamp.Layout.Row = 1;
            app.apd1LoadedLamp.Layout.Column = 7;
            app.apd1LoadedLamp.Color = [1 0 0];

            % Create dataFielesListBox
            app.dataFielesListBox = uilistbox(app.GridLayout);
            app.dataFielesListBox.Items = {'No Data file selected!'};
            app.dataFielesListBox.Multiselect = 'on';
            app.dataFielesListBox.ValueChangedFcn = createCallbackFcn(app, @dataFielesListBoxValueChanged, true);
            app.dataFielesListBox.Layout.Row = [2 8];
            app.dataFielesListBox.Layout.Column = [1 3];
            app.dataFielesListBox.Value = {'No Data file selected!'};

            % Create dataFilesLabel
            app.dataFilesLabel = uilabel(app.GridLayout);
            app.dataFilesLabel.HorizontalAlignment = 'right';
            app.dataFilesLabel.Layout.Row = 1;
            app.dataFilesLabel.Layout.Column = 1;
            app.dataFilesLabel.Text = 'Data Files';

            % Create DeConSettingTab
            app.DeConSettingTab = uitab(app.TabGroup);
            app.DeConSettingTab.Title = 'DeCon Setting';

            % Create GridLayout2
            app.GridLayout2 = uigridlayout(app.DeConSettingTab);
            app.GridLayout2.ColumnWidth = {'0.5x', '0.5x', '0.5x', '0.5x', '1x', '1x', '1x', '1x', '1x', '1x'};
            app.GridLayout2.RowHeight = {'1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1.15x'};
            app.GridLayout2.ColumnSpacing = 3.5;
            app.GridLayout2.RowSpacing = 2;
            app.GridLayout2.Padding = [3.5 2 3.5 2];
            app.GridLayout2.BusyAction = 'cancel';

            % Create UIAxes_ch1DataFig
            app.UIAxes_ch1DataFig = uiaxes(app.GridLayout2);
            title(app.UIAxes_ch1DataFig, 'Ch1')
            xlabel(app.UIAxes_ch1DataFig, 'Points')
            ylabel(app.UIAxes_ch1DataFig, 'Voltage')
            app.UIAxes_ch1DataFig.Box = 'on';
            app.UIAxes_ch1DataFig.XGrid = 'on';
            app.UIAxes_ch1DataFig.YGrid = 'on';
            app.UIAxes_ch1DataFig.Layout.Row = [2 9];
            app.UIAxes_ch1DataFig.Layout.Column = [5 7];

            % Create UIAxes_ch2DataFig
            app.UIAxes_ch2DataFig = uiaxes(app.GridLayout2);
            title(app.UIAxes_ch2DataFig, 'Ch2')
            xlabel(app.UIAxes_ch2DataFig, 'Points')
            ylabel(app.UIAxes_ch2DataFig, 'Voltage')
            app.UIAxes_ch2DataFig.Box = 'on';
            app.UIAxes_ch2DataFig.XGrid = 'on';
            app.UIAxes_ch2DataFig.YGrid = 'on';
            app.UIAxes_ch2DataFig.Layout.Row = [10 17];
            app.UIAxes_ch2DataFig.Layout.Column = [5 7];

            % Create UIAxes_ch1BgRemoveFig
            app.UIAxes_ch1BgRemoveFig = uiaxes(app.GridLayout2);
            title(app.UIAxes_ch1BgRemoveFig, 'background removed ch1')
            xlabel(app.UIAxes_ch1BgRemoveFig, 'Points')
            ylabel(app.UIAxes_ch1BgRemoveFig, 'Voltage')
            app.UIAxes_ch1BgRemoveFig.Box = 'on';
            app.UIAxes_ch1BgRemoveFig.XGrid = 'on';
            app.UIAxes_ch1BgRemoveFig.YGrid = 'on';
            app.UIAxes_ch1BgRemoveFig.Layout.Row = [2 9];
            app.UIAxes_ch1BgRemoveFig.Layout.Column = [8 10];

            % Create UIAxes_ch3DataFig
            app.UIAxes_ch3DataFig = uiaxes(app.GridLayout2);
            title(app.UIAxes_ch3DataFig, 'Ch3')
            xlabel(app.UIAxes_ch3DataFig, 'Points')
            ylabel(app.UIAxes_ch3DataFig, 'Voltage')
            app.UIAxes_ch3DataFig.Box = 'on';
            app.UIAxes_ch3DataFig.XGrid = 'on';
            app.UIAxes_ch3DataFig.YGrid = 'on';
            app.UIAxes_ch3DataFig.Layout.Row = [18 25];
            app.UIAxes_ch3DataFig.Layout.Column = [5 7];

            % Create UIAxes_ch2BgRemoveFig
            app.UIAxes_ch2BgRemoveFig = uiaxes(app.GridLayout2);
            title(app.UIAxes_ch2BgRemoveFig, 'background removed ch2')
            xlabel(app.UIAxes_ch2BgRemoveFig, 'Points')
            ylabel(app.UIAxes_ch2BgRemoveFig, 'Voltage')
            app.UIAxes_ch2BgRemoveFig.Box = 'on';
            app.UIAxes_ch2BgRemoveFig.XGrid = 'on';
            app.UIAxes_ch2BgRemoveFig.YGrid = 'on';
            app.UIAxes_ch2BgRemoveFig.Layout.Row = [10 17];
            app.UIAxes_ch2BgRemoveFig.Layout.Column = [8 10];

            % Create UIAxes_ch3BgRemoveFig
            app.UIAxes_ch3BgRemoveFig = uiaxes(app.GridLayout2);
            title(app.UIAxes_ch3BgRemoveFig, 'background removed ch3')
            xlabel(app.UIAxes_ch3BgRemoveFig, 'Points')
            ylabel(app.UIAxes_ch3BgRemoveFig, 'Voltage')
            app.UIAxes_ch3BgRemoveFig.Box = 'on';
            app.UIAxes_ch3BgRemoveFig.XGrid = 'on';
            app.UIAxes_ch3BgRemoveFig.YGrid = 'on';
            app.UIAxes_ch3BgRemoveFig.Layout.Row = [18 25];
            app.UIAxes_ch3BgRemoveFig.Layout.Column = [8 10];

            % Create RemoveSaturationButton
            app.RemoveSaturationButton = uibutton(app.GridLayout2, 'push');
            app.RemoveSaturationButton.ButtonPushedFcn = createCallbackFcn(app, @RemoveSaturationButtonPushed, true);
            app.RemoveSaturationButton.FontSize = 15;
            app.RemoveSaturationButton.FontWeight = 'bold';
            app.RemoveSaturationButton.Enable = 'off';
            app.RemoveSaturationButton.Layout.Row = 26;
            app.RemoveSaturationButton.Layout.Column = 6;
            app.RemoveSaturationButton.Text = 'Remove Saturation';

            % Create TruncateButton
            app.TruncateButton = uibutton(app.GridLayout2, 'push');
            app.TruncateButton.ButtonPushedFcn = createCallbackFcn(app, @TruncateButtonPushed, true);
            app.TruncateButton.FontSize = 18;
            app.TruncateButton.FontWeight = 'bold';
            app.TruncateButton.Enable = 'off';
            app.TruncateButton.Layout.Row = 26;
            app.TruncateButton.Layout.Column = 9;
            app.TruncateButton.Text = 'Truncate';

            % Create AverageDataButton
            app.AverageDataButton = uibutton(app.GridLayout2, 'push');
            app.AverageDataButton.ButtonPushedFcn = createCallbackFcn(app, @AverageDataButtonPushed, true);
            app.AverageDataButton.FontSize = 18;
            app.AverageDataButton.FontWeight = 'bold';
            app.AverageDataButton.Enable = 'off';
            app.AverageDataButton.Layout.Row = 26;
            app.AverageDataButton.Layout.Column = 7;
            app.AverageDataButton.Text = 'Average Data';

            % Create DeconButton
            app.DeconButton = uibutton(app.GridLayout2, 'push');
            app.DeconButton.ButtonPushedFcn = createCallbackFcn(app, @DeconButtonPushed, true);
            app.DeconButton.FontSize = 18;
            app.DeconButton.FontWeight = 'bold';
            app.DeconButton.Enable = 'off';
            app.DeconButton.Layout.Row = 26;
            app.DeconButton.Layout.Column = 10;
            app.DeconButton.Text = 'Decon';

            % Create nmlaserreprateEditFieldLabel
            app.nmlaserreprateEditFieldLabel = uilabel(app.GridLayout2);
            app.nmlaserreprateEditFieldLabel.HorizontalAlignment = 'center';
            app.nmlaserreprateEditFieldLabel.FontSize = 14;
            app.nmlaserreprateEditFieldLabel.FontWeight = 'bold';
            app.nmlaserreprateEditFieldLabel.FontColor = [0 0 1];
            app.nmlaserreprateEditFieldLabel.Layout.Row = 1;
            app.nmlaserreprateEditFieldLabel.Layout.Column = [1 2];
            app.nmlaserreprateEditFieldLabel.Text = '355nm laser rep rate';

            % Create laserRepRateEditField
            app.laserRepRateEditField = uieditfield(app.GridLayout2, 'numeric');
            app.laserRepRateEditField.Limits = [0 Inf];
            app.laserRepRateEditField.FontSize = 14;
            app.laserRepRateEditField.FontWeight = 'bold';
            app.laserRepRateEditField.Layout.Row = 1;
            app.laserRepRateEditField.Layout.Column = [3 4];

            % Create UpsampleButton
            app.UpsampleButton = uibutton(app.GridLayout2, 'push');
            app.UpsampleButton.ButtonPushedFcn = createCallbackFcn(app, @UpsampleButtonPushed, true);
            app.UpsampleButton.FontSize = 18;
            app.UpsampleButton.FontWeight = 'bold';
            app.UpsampleButton.Enable = 'off';
            app.UpsampleButton.Layout.Row = 26;
            app.UpsampleButton.Layout.Column = 5;
            app.UpsampleButton.Text = 'Upsample';

            % Create RemoveBGButton
            app.RemoveBGButton = uibutton(app.GridLayout2, 'push');
            app.RemoveBGButton.ButtonPushedFcn = createCallbackFcn(app, @RemoveBGButtonPushed, true);
            app.RemoveBGButton.FontSize = 18;
            app.RemoveBGButton.FontWeight = 'bold';
            app.RemoveBGButton.Enable = 'off';
            app.RemoveBGButton.Layout.Row = 26;
            app.RemoveBGButton.Layout.Column = 8;
            app.RemoveBGButton.Text = 'Remove BG';

            % Create LoadDataButton
            app.LoadDataButton = uibutton(app.GridLayout2, 'push');
            app.LoadDataButton.ButtonPushedFcn = createCallbackFcn(app, @LoadDataButtonPushed, true);
            app.LoadDataButton.FontSize = 18;
            app.LoadDataButton.FontWeight = 'bold';
            app.LoadDataButton.Layout.Row = 26;
            app.LoadDataButton.Layout.Column = [3 4];
            app.LoadDataButton.Text = 'Load Data';

            % Create channelWidthEditField
            app.channelWidthEditField = uieditfield(app.GridLayout2, 'numeric');
            app.channelWidthEditField.Limits = [0 Inf];
            app.channelWidthEditField.ValueDisplayFormat = '%3.3f';
            app.channelWidthEditField.ValueChangedFcn = createCallbackFcn(app, @channelWidthEditFieldValueChanged, true);
            app.channelWidthEditField.Enable = 'off';
            app.channelWidthEditField.Layout.Row = 14;
            app.channelWidthEditField.Layout.Column = [3 4];
            app.channelWidthEditField.Value = 154;

            % Create ChannelWidthnsLabel
            app.ChannelWidthnsLabel = uilabel(app.GridLayout2);
            app.ChannelWidthnsLabel.HorizontalAlignment = 'right';
            app.ChannelWidthnsLabel.Enable = 'off';
            app.ChannelWidthnsLabel.Layout.Row = 14;
            app.ChannelWidthnsLabel.Layout.Column = [1 2];
            app.ChannelWidthnsLabel.Text = 'Channel Width (ns)';

            % Create LaguerreOrderEditField
            app.LaguerreOrderEditField = uieditfield(app.GridLayout2, 'numeric');
            app.LaguerreOrderEditField.Limits = [1 30];
            app.LaguerreOrderEditField.ValueDisplayFormat = '%.0f';
            app.LaguerreOrderEditField.Enable = 'off';
            app.LaguerreOrderEditField.Layout.Row = 18;
            app.LaguerreOrderEditField.Layout.Column = [3 4];
            app.LaguerreOrderEditField.Value = 12;

            % Create LaguerreOrderEditFieldLabel
            app.LaguerreOrderEditFieldLabel = uilabel(app.GridLayout2);
            app.LaguerreOrderEditFieldLabel.HorizontalAlignment = 'right';
            app.LaguerreOrderEditFieldLabel.Enable = 'off';
            app.LaguerreOrderEditFieldLabel.Layout.Row = 18;
            app.LaguerreOrderEditFieldLabel.Layout.Column = [1 2];
            app.LaguerreOrderEditFieldLabel.Text = 'Laguerre Order';

            % Create dataIgnoreLowEditField
            app.dataIgnoreLowEditField = uieditfield(app.GridLayout2, 'numeric');
            app.dataIgnoreLowEditField.RoundFractionalValues = 'on';
            app.dataIgnoreLowEditField.ValueChangedFcn = createCallbackFcn(app, @dataIgnoreLowEditFieldValueChanged, true);
            app.dataIgnoreLowEditField.Enable = 'off';
            app.dataIgnoreLowEditField.Layout.Row = 16;
            app.dataIgnoreLowEditField.Layout.Column = 2;
            app.dataIgnoreLowEditField.Value = 580;

            % Create DcLowEditFieldLabel
            app.DcLowEditFieldLabel = uilabel(app.GridLayout2);
            app.DcLowEditFieldLabel.HorizontalAlignment = 'right';
            app.DcLowEditFieldLabel.Enable = 'off';
            app.DcLowEditFieldLabel.Layout.Row = 16;
            app.DcLowEditFieldLabel.Layout.Column = 1;
            app.DcLowEditFieldLabel.Text = 'Low';

            % Create dataIgnoreHighEditField
            app.dataIgnoreHighEditField = uieditfield(app.GridLayout2, 'numeric');
            app.dataIgnoreHighEditField.RoundFractionalValues = 'on';
            app.dataIgnoreHighEditField.ValueChangedFcn = createCallbackFcn(app, @dataIgnoreHighEditFieldValueChanged, true);
            app.dataIgnoreHighEditField.Enable = 'off';
            app.dataIgnoreHighEditField.Layout.Row = 16;
            app.dataIgnoreHighEditField.Layout.Column = 4;
            app.dataIgnoreHighEditField.Value = 660;

            % Create DcHighEditFieldLabel
            app.DcHighEditFieldLabel = uilabel(app.GridLayout2);
            app.DcHighEditFieldLabel.HorizontalAlignment = 'right';
            app.DcHighEditFieldLabel.Enable = 'off';
            app.DcHighEditFieldLabel.Layout.Row = 16;
            app.DcHighEditFieldLabel.Layout.Column = 3;
            app.DcHighEditFieldLabel.Text = 'High';

            % Create BgLowEditField
            app.BgLowEditField = uieditfield(app.GridLayout2, 'numeric');
            app.BgLowEditField.ValueChangedFcn = createCallbackFcn(app, @BgLowEditFieldValueChanged, true);
            app.BgLowEditField.Enable = 'off';
            app.BgLowEditField.Layout.Row = 11;
            app.BgLowEditField.Layout.Column = 2;
            app.BgLowEditField.Value = 100;

            % Create BgLowEditFieldLabel
            app.BgLowEditFieldLabel = uilabel(app.GridLayout2);
            app.BgLowEditFieldLabel.HorizontalAlignment = 'right';
            app.BgLowEditFieldLabel.Enable = 'off';
            app.BgLowEditFieldLabel.Layout.Row = 11;
            app.BgLowEditFieldLabel.Layout.Column = 1;
            app.BgLowEditFieldLabel.Text = 'Low';

            % Create BgHighEditField
            app.BgHighEditField = uieditfield(app.GridLayout2, 'numeric');
            app.BgHighEditField.ValueChangedFcn = createCallbackFcn(app, @BgHighEditFieldValueChanged, true);
            app.BgHighEditField.Enable = 'off';
            app.BgHighEditField.Layout.Row = 11;
            app.BgHighEditField.Layout.Column = 4;
            app.BgHighEditField.Value = 175;

            % Create BgHighEditFieldLabel
            app.BgHighEditFieldLabel = uilabel(app.GridLayout2);
            app.BgHighEditFieldLabel.HorizontalAlignment = 'right';
            app.BgHighEditFieldLabel.Enable = 'off';
            app.BgHighEditFieldLabel.Layout.Row = 11;
            app.BgHighEditFieldLabel.Layout.Column = 3;
            app.BgHighEditFieldLabel.Text = 'High';

            % Create BgThreholdEditField
            app.BgThreholdEditField = uieditfield(app.GridLayout2, 'numeric');
            app.BgThreholdEditField.Enable = 'off';
            app.BgThreholdEditField.Layout.Row = 12;
            app.BgThreholdEditField.Layout.Column = 4;
            app.BgThreholdEditField.Value = 0.003;

            % Create BgThreholdEditFieldLabel
            app.BgThreholdEditFieldLabel = uilabel(app.GridLayout2);
            app.BgThreholdEditFieldLabel.HorizontalAlignment = 'right';
            app.BgThreholdEditFieldLabel.Enable = 'off';
            app.BgThreholdEditFieldLabel.Layout.Row = 12;
            app.BgThreholdEditFieldLabel.Layout.Column = 3;
            app.BgThreholdEditFieldLabel.Text = 'Threhold';

            % Create TimeResolutionnsEditField
            app.TimeResolutionnsEditField = uieditfield(app.GridLayout2, 'numeric');
            app.TimeResolutionnsEditField.Enable = 'off';
            app.TimeResolutionnsEditField.Layout.Row = 7;
            app.TimeResolutionnsEditField.Layout.Column = [3 4];
            app.TimeResolutionnsEditField.Value = 0.4;

            % Create TimeResolutionnsEditFieldLabel
            app.TimeResolutionnsEditFieldLabel = uilabel(app.GridLayout2);
            app.TimeResolutionnsEditFieldLabel.HorizontalAlignment = 'right';
            app.TimeResolutionnsEditFieldLabel.Enable = 'off';
            app.TimeResolutionnsEditFieldLabel.Layout.Row = 7;
            app.TimeResolutionnsEditFieldLabel.Layout.Column = [1 2];
            app.TimeResolutionnsEditFieldLabel.Text = 'Time Resolution (ns)';

            % Create DataAverageDropDownLabel
            app.DataAverageDropDownLabel = uilabel(app.GridLayout2);
            app.DataAverageDropDownLabel.HorizontalAlignment = 'right';
            app.DataAverageDropDownLabel.Enable = 'off';
            app.DataAverageDropDownLabel.Layout.Row = 5;
            app.DataAverageDropDownLabel.Layout.Column = [1 2];
            app.DataAverageDropDownLabel.Text = 'Data Average';

            % Create DataAverageDropDown
            app.DataAverageDropDown = uidropdown(app.GridLayout2);
            app.DataAverageDropDown.Items = {'1'};
            app.DataAverageDropDown.Enable = 'off';
            app.DataAverageDropDown.Layout.Row = 5;
            app.DataAverageDropDown.Layout.Column = [3 4];
            app.DataAverageDropDown.Value = '1';

            % Create LowThresholdEditField
            app.LowThresholdEditField = uieditfield(app.GridLayout2, 'numeric');
            app.LowThresholdEditField.Enable = 'off';
            app.LowThresholdEditField.Layout.Row = 3;
            app.LowThresholdEditField.Layout.Column = 2;
            app.LowThresholdEditField.Value = 0.03;

            % Create LowThresholdEditFieldLabel
            app.LowThresholdEditFieldLabel = uilabel(app.GridLayout2);
            app.LowThresholdEditFieldLabel.HorizontalAlignment = 'right';
            app.LowThresholdEditFieldLabel.Enable = 'off';
            app.LowThresholdEditFieldLabel.Layout.Row = 3;
            app.LowThresholdEditFieldLabel.Layout.Column = 1;
            app.LowThresholdEditFieldLabel.Text = 'Low Threshold';

            % Create HighThresholdEditField
            app.HighThresholdEditField = uieditfield(app.GridLayout2, 'numeric');
            app.HighThresholdEditField.Limits = [0 1.7];
            app.HighThresholdEditField.Enable = 'off';
            app.HighThresholdEditField.Layout.Row = 3;
            app.HighThresholdEditField.Layout.Column = 4;
            app.HighThresholdEditField.Value = 1.4;

            % Create HighThresholdEditFieldLabel
            app.HighThresholdEditFieldLabel = uilabel(app.GridLayout2);
            app.HighThresholdEditFieldLabel.HorizontalAlignment = 'right';
            app.HighThresholdEditFieldLabel.Enable = 'off';
            app.HighThresholdEditFieldLabel.Layout.Row = 3;
            app.HighThresholdEditFieldLabel.Layout.Column = 3;
            app.HighThresholdEditFieldLabel.Text = 'HighThreshold';

            % Create SignalConditioningPanel
            app.SignalConditioningPanel = uipanel(app.GridLayout2);
            app.SignalConditioningPanel.ForegroundColor = [1 0 0];
            app.SignalConditioningPanel.TitlePosition = 'centertop';
            app.SignalConditioningPanel.Title = 'Signal Conditioning';
            app.SignalConditioningPanel.Layout.Row = 2;
            app.SignalConditioningPanel.Layout.Column = [1 4];
            app.SignalConditioningPanel.FontWeight = 'bold';
            app.SignalConditioningPanel.FontSize = 14;

            % Create DataAveragingPanel
            app.DataAveragingPanel = uipanel(app.GridLayout2);
            app.DataAveragingPanel.ForegroundColor = [1 0 0];
            app.DataAveragingPanel.TitlePosition = 'centertop';
            app.DataAveragingPanel.Title = 'Data Averaging';
            app.DataAveragingPanel.Layout.Row = 4;
            app.DataAveragingPanel.Layout.Column = [1 4];
            app.DataAveragingPanel.FontWeight = 'bold';
            app.DataAveragingPanel.FontSize = 14;

            % Create DigitizerSettingPanel
            app.DigitizerSettingPanel = uipanel(app.GridLayout2);
            app.DigitizerSettingPanel.ForegroundColor = [1 0 0];
            app.DigitizerSettingPanel.TitlePosition = 'centertop';
            app.DigitizerSettingPanel.Title = 'Digitizer Setting';
            app.DigitizerSettingPanel.Layout.Row = 6;
            app.DigitizerSettingPanel.Layout.Column = [1 4];
            app.DigitizerSettingPanel.FontWeight = 'bold';
            app.DigitizerSettingPanel.FontSize = 14;

            % Create channel1CheckBox
            app.channel1CheckBox = uicheckbox(app.GridLayout2);
            app.channel1CheckBox.Enable = 'off';
            app.channel1CheckBox.Text = 'Channel 1';
            app.channel1CheckBox.Layout.Row = 9;
            app.channel1CheckBox.Layout.Column = 1;
            app.channel1CheckBox.Value = true;

            % Create channel2CheckBox
            app.channel2CheckBox = uicheckbox(app.GridLayout2);
            app.channel2CheckBox.Enable = 'off';
            app.channel2CheckBox.Text = 'Channel 2';
            app.channel2CheckBox.Layout.Row = 9;
            app.channel2CheckBox.Layout.Column = 2;
            app.channel2CheckBox.Value = true;

            % Create channel3CheckBox
            app.channel3CheckBox = uicheckbox(app.GridLayout2);
            app.channel3CheckBox.Enable = 'off';
            app.channel3CheckBox.Text = 'Channel 3';
            app.channel3CheckBox.Layout.Row = 9;
            app.channel3CheckBox.Layout.Column = 3;
            app.channel3CheckBox.Value = true;

            % Create ChannelSelectionPanel
            app.ChannelSelectionPanel = uipanel(app.GridLayout2);
            app.ChannelSelectionPanel.ForegroundColor = [1 0 0];
            app.ChannelSelectionPanel.TitlePosition = 'centertop';
            app.ChannelSelectionPanel.Title = 'Channel Selection';
            app.ChannelSelectionPanel.Layout.Row = 8;
            app.ChannelSelectionPanel.Layout.Column = [1 4];
            app.ChannelSelectionPanel.FontWeight = 'bold';
            app.ChannelSelectionPanel.FontSize = 14;

            % Create BackgroudRemovalPanel
            app.BackgroudRemovalPanel = uipanel(app.GridLayout2);
            app.BackgroudRemovalPanel.ForegroundColor = [1 0 0];
            app.BackgroudRemovalPanel.TitlePosition = 'centertop';
            app.BackgroudRemovalPanel.Title = 'Backgroud Removal';
            app.BackgroudRemovalPanel.Layout.Row = 10;
            app.BackgroudRemovalPanel.Layout.Column = [1 4];
            app.BackgroudRemovalPanel.FontWeight = 'bold';
            app.BackgroudRemovalPanel.FontSize = 14;

            % Create IgnoreDataRangePanel
            app.IgnoreDataRangePanel = uipanel(app.GridLayout2);
            app.IgnoreDataRangePanel.ForegroundColor = [1 0 0];
            app.IgnoreDataRangePanel.TitlePosition = 'centertop';
            app.IgnoreDataRangePanel.Title = 'Ignore Data Range';
            app.IgnoreDataRangePanel.Layout.Row = 15;
            app.IgnoreDataRangePanel.Layout.Column = [1 4];
            app.IgnoreDataRangePanel.FontWeight = 'bold';
            app.IgnoreDataRangePanel.FontSize = 14;

            % Create MultiExponentialSettingPanel
            app.MultiExponentialSettingPanel = uipanel(app.GridLayout2);
            app.MultiExponentialSettingPanel.ForegroundColor = [1 0 0];
            app.MultiExponentialSettingPanel.TitlePosition = 'centertop';
            app.MultiExponentialSettingPanel.Title = 'Multi-Exponential Setting';
            app.MultiExponentialSettingPanel.Layout.Row = 23;
            app.MultiExponentialSettingPanel.Layout.Column = [1 4];
            app.MultiExponentialSettingPanel.FontWeight = 'bold';
            app.MultiExponentialSettingPanel.FontSize = 14;

            % Create ofexponentialsDropDownLabel
            app.ofexponentialsDropDownLabel = uilabel(app.GridLayout2);
            app.ofexponentialsDropDownLabel.HorizontalAlignment = 'right';
            app.ofexponentialsDropDownLabel.Enable = 'off';
            app.ofexponentialsDropDownLabel.Layout.Row = 24;
            app.ofexponentialsDropDownLabel.Layout.Column = [2 3];
            app.ofexponentialsDropDownLabel.Text = '# of exponentials';

            % Create exponentialsDropDown
            app.exponentialsDropDown = uidropdown(app.GridLayout2);
            app.exponentialsDropDown.Items = {'1', '2', '3', '4'};
            app.exponentialsDropDown.Enable = 'off';
            app.exponentialsDropDown.Layout.Row = 24;
            app.exponentialsDropDown.Layout.Column = 4;
            app.exponentialsDropDown.Value = '3';

            % Create LaguerreSettingPanel
            app.LaguerreSettingPanel = uipanel(app.GridLayout2);
            app.LaguerreSettingPanel.ForegroundColor = [1 0 0];
            app.LaguerreSettingPanel.TitlePosition = 'centertop';
            app.LaguerreSettingPanel.Title = 'Laguerre Setting';
            app.LaguerreSettingPanel.Layout.Row = 17;
            app.LaguerreSettingPanel.Layout.Column = [1 4];
            app.LaguerreSettingPanel.FontWeight = 'bold';
            app.LaguerreSettingPanel.FontSize = 14;

            % Create Ch1AlphaValueDropDownLabel
            app.Ch1AlphaValueDropDownLabel = uilabel(app.GridLayout2);
            app.Ch1AlphaValueDropDownLabel.HorizontalAlignment = 'right';
            app.Ch1AlphaValueDropDownLabel.Enable = 'off';
            app.Ch1AlphaValueDropDownLabel.Layout.Row = 19;
            app.Ch1AlphaValueDropDownLabel.Layout.Column = [1 2];
            app.Ch1AlphaValueDropDownLabel.Text = 'Ch1 Alpha Value';

            % Create Ch1AlphaValueDropDown
            app.Ch1AlphaValueDropDown = uidropdown(app.GridLayout2);
            app.Ch1AlphaValueDropDown.Items = {'0.880 (0.5-5 ns)', '0.916 (0.7-8 ns)', '0.930 (1-8 ns)', '0.965 (2-16 ns)'};
            app.Ch1AlphaValueDropDown.ItemsData = {'0.880', '0.916', '0.930', '0.965', ''};
            app.Ch1AlphaValueDropDown.Enable = 'off';
            app.Ch1AlphaValueDropDown.Layout.Row = 19;
            app.Ch1AlphaValueDropDown.Layout.Column = [3 4];
            app.Ch1AlphaValueDropDown.Value = '0.880';

            % Create Ch2AlphaValueDropDownLabel
            app.Ch2AlphaValueDropDownLabel = uilabel(app.GridLayout2);
            app.Ch2AlphaValueDropDownLabel.HorizontalAlignment = 'right';
            app.Ch2AlphaValueDropDownLabel.Enable = 'off';
            app.Ch2AlphaValueDropDownLabel.Layout.Row = 20;
            app.Ch2AlphaValueDropDownLabel.Layout.Column = [1 2];
            app.Ch2AlphaValueDropDownLabel.Text = 'Ch2 Alpha Value';

            % Create Ch2AlphaValueDropDown
            app.Ch2AlphaValueDropDown = uidropdown(app.GridLayout2);
            app.Ch2AlphaValueDropDown.Items = {'0.880 (0.5-5 ns)', '0.916 (0.7-8 ns)', '0.930 (1-8 ns)', '0.965 (2-16 ns)'};
            app.Ch2AlphaValueDropDown.ItemsData = {'0.880', '0.916', '0.930', '0.965', ''};
            app.Ch2AlphaValueDropDown.Enable = 'off';
            app.Ch2AlphaValueDropDown.Layout.Row = 20;
            app.Ch2AlphaValueDropDown.Layout.Column = [3 4];
            app.Ch2AlphaValueDropDown.Value = '0.916';

            % Create Ch3AlphaValueDropDownLabel
            app.Ch3AlphaValueDropDownLabel = uilabel(app.GridLayout2);
            app.Ch3AlphaValueDropDownLabel.HorizontalAlignment = 'right';
            app.Ch3AlphaValueDropDownLabel.Enable = 'off';
            app.Ch3AlphaValueDropDownLabel.Layout.Row = 21;
            app.Ch3AlphaValueDropDownLabel.Layout.Column = [1 2];
            app.Ch3AlphaValueDropDownLabel.Text = 'Ch3 Alpha Value';

            % Create Ch3AlphaValueDropDown
            app.Ch3AlphaValueDropDown = uidropdown(app.GridLayout2);
            app.Ch3AlphaValueDropDown.Items = {'0.880 (0.5-5 ns)', '0.916 (0.7-8 ns)', '0.930 (1-8 ns)', '0.965 (2-16 ns)'};
            app.Ch3AlphaValueDropDown.ItemsData = {'0.880', '0.916', '0.930', '0.965', ''};
            app.Ch3AlphaValueDropDown.Enable = 'off';
            app.Ch3AlphaValueDropDown.Layout.Row = 21;
            app.Ch3AlphaValueDropDown.Layout.Column = [3 4];
            app.Ch3AlphaValueDropDown.Value = '0.916';

            % Create TruncationPanel
            app.TruncationPanel = uipanel(app.GridLayout2);
            app.TruncationPanel.ForegroundColor = [1 0 0];
            app.TruncationPanel.TitlePosition = 'centertop';
            app.TruncationPanel.Title = 'Truncation';
            app.TruncationPanel.Layout.Row = 13;
            app.TruncationPanel.Layout.Column = [1 4];
            app.TruncationPanel.FontWeight = 'bold';
            app.TruncationPanel.FontSize = 14;

            % Create LogScaleCheckBox
            app.LogScaleCheckBox = uicheckbox(app.GridLayout2);
            app.LogScaleCheckBox.Text = 'Log Scale';
            app.LogScaleCheckBox.FontWeight = 'bold';
            app.LogScaleCheckBox.FontColor = [0 0 1];
            app.LogScaleCheckBox.Layout.Row = 1;
            app.LogScaleCheckBox.Layout.Column = 10;

            % Create RunEXPCheckBox
            app.RunEXPCheckBox = uicheckbox(app.GridLayout2);
            app.RunEXPCheckBox.ValueChangedFcn = createCallbackFcn(app, @RunEXPCheckBoxValueChanged, true);
            app.RunEXPCheckBox.Text = 'Run EXP';
            app.RunEXPCheckBox.Layout.Row = 24;
            app.RunEXPCheckBox.Layout.Column = 1;

            % Create FittingResultTab
            app.FittingResultTab = uitab(app.TabGroup);
            app.FittingResultTab.Title = 'Fitting Result';

            % Create GridLayout3
            app.GridLayout3 = uigridlayout(app.FittingResultTab);
            app.GridLayout3.ColumnWidth = {'1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x'};
            app.GridLayout3.RowHeight = {22, '60x', '25x', '15x', '10x'};
            app.GridLayout3.RowSpacing = 2.75;
            app.GridLayout3.Padding = [10 2.75 10 2.75];

            % Create UIAxesFitting
            app.UIAxesFitting = uiaxes(app.GridLayout3);
            title(app.UIAxesFitting, 'Fitting')
            xlabel(app.UIAxesFitting, 'Time(ns)')
            ylabel(app.UIAxesFitting, 'Y')
            app.UIAxesFitting.XTickLabelRotation = 0;
            app.UIAxesFitting.YTickLabelRotation = 0;
            app.UIAxesFitting.ZTickLabelRotation = 0;
            app.UIAxesFitting.Box = 'on';
            app.UIAxesFitting.XGrid = 'on';
            app.UIAxesFitting.YGrid = 'on';
            app.UIAxesFitting.Layout.Row = 2;
            app.UIAxesFitting.Layout.Column = [1 6];

            % Create UIAxesAutoCo
            app.UIAxesAutoCo = uiaxes(app.GridLayout3);
            title(app.UIAxesAutoCo, 'Auto Correlation')
            ylabel(app.UIAxesAutoCo, 'Y')
            app.UIAxesAutoCo.XTickLabelRotation = 0;
            app.UIAxesAutoCo.YTickLabelRotation = 0;
            app.UIAxesAutoCo.ZTickLabelRotation = 0;
            app.UIAxesAutoCo.Box = 'on';
            app.UIAxesAutoCo.XGrid = 'on';
            app.UIAxesAutoCo.YGrid = 'on';
            app.UIAxesAutoCo.Layout.Row = [4 5];
            app.UIAxesAutoCo.Layout.Column = [1 6];

            % Create UIAxesResidue
            app.UIAxesResidue = uiaxes(app.GridLayout3);
            title(app.UIAxesResidue, 'Residue')
            ylabel(app.UIAxesResidue, 'Y')
            app.UIAxesResidue.XTickLabelRotation = 0;
            app.UIAxesResidue.YTickLabelRotation = 0;
            app.UIAxesResidue.ZTickLabelRotation = 0;
            app.UIAxesResidue.Box = 'on';
            app.UIAxesResidue.XGrid = 'on';
            app.UIAxesResidue.YGrid = 'on';
            app.UIAxesResidue.Layout.Row = 3;
            app.UIAxesResidue.Layout.Column = [1 6];

            % Create UIAxesDeconResult
            app.UIAxesDeconResult = uiaxes(app.GridLayout3);
            title(app.UIAxesDeconResult, 'Utility Plot')
            xlabel(app.UIAxesDeconResult, 'Time(ns)')
            ylabel(app.UIAxesDeconResult, 'Y')
            app.UIAxesDeconResult.Box = 'on';
            app.UIAxesDeconResult.XGrid = 'on';
            app.UIAxesDeconResult.YGrid = 'on';
            app.UIAxesDeconResult.Layout.Row = [2 3];
            app.UIAxesDeconResult.Layout.Column = [7 12];

            % Create channelDropDown
            app.channelDropDown = uidropdown(app.GridLayout3);
            app.channelDropDown.ValueChangedFcn = createCallbackFcn(app, @channelDropDownValueChanged, true);
            app.channelDropDown.Layout.Row = 1;
            app.channelDropDown.Layout.Column = 2;

            % Create channelLabel
            app.channelLabel = uilabel(app.GridLayout3);
            app.channelLabel.HorizontalAlignment = 'right';
            app.channelLabel.FontWeight = 'bold';
            app.channelLabel.Layout.Row = 1;
            app.channelLabel.Layout.Column = 1;
            app.channelLabel.Text = 'Channel #';

            % Create PlotDropDown
            app.PlotDropDown = uidropdown(app.GridLayout3);
            app.PlotDropDown.Items = {'Lifetime', 'Laguerre basis', 'Raw Intensity', 'Decon Intensity', 'Gain', 'Raw Waveforms', 'Residual Max', 'SNR', 'LG vs mExp', 'Square Error', 'Shift'};
            app.PlotDropDown.ValueChangedFcn = createCallbackFcn(app, @PlotDropDownValueChanged, true);
            app.PlotDropDown.BackgroundColor = [1 1 1];
            app.PlotDropDown.Layout.Row = 1;
            app.PlotDropDown.Layout.Column = [8 9];
            app.PlotDropDown.Value = 'Lifetime';

            % Create xMinEditField
            app.xMinEditField = uieditfield(app.GridLayout3, 'numeric');
            app.xMinEditField.ValueChangedFcn = createCallbackFcn(app, @xMinEditFieldValueChanged, true);
            app.xMinEditField.Layout.Row = 1;
            app.xMinEditField.Layout.Column = 4;

            % Create xMinEditFieldLabel
            app.xMinEditFieldLabel = uilabel(app.GridLayout3);
            app.xMinEditFieldLabel.HorizontalAlignment = 'right';
            app.xMinEditFieldLabel.Layout.Row = 1;
            app.xMinEditFieldLabel.Layout.Column = 3;
            app.xMinEditFieldLabel.Text = 'xMin';

            % Create xMaxEditField
            app.xMaxEditField = uieditfield(app.GridLayout3, 'numeric');
            app.xMaxEditField.ValueDisplayFormat = '%5.2f';
            app.xMaxEditField.ValueChangedFcn = createCallbackFcn(app, @xMinEditFieldValueChanged, true);
            app.xMaxEditField.Layout.Row = 1;
            app.xMaxEditField.Layout.Column = 6;
            app.xMaxEditField.Value = 154;

            % Create xMaxEditFieldLabel
            app.xMaxEditFieldLabel = uilabel(app.GridLayout3);
            app.xMaxEditFieldLabel.HorizontalAlignment = 'right';
            app.xMaxEditFieldLabel.Layout.Row = 1;
            app.xMaxEditFieldLabel.Layout.Column = 5;
            app.xMaxEditFieldLabel.Text = 'xMax';

            % Create SaveDataButton
            app.SaveDataButton = uibutton(app.GridLayout3, 'push');
            app.SaveDataButton.ButtonPushedFcn = createCallbackFcn(app, @SaveDataButtonPushed, true);
            app.SaveDataButton.Layout.Row = 5;
            app.SaveDataButton.Layout.Column = [11 12];
            app.SaveDataButton.Text = 'Save Data';

            % Create PlotLabel
            app.PlotLabel = uilabel(app.GridLayout3);
            app.PlotLabel.HorizontalAlignment = 'right';
            app.PlotLabel.FontWeight = 'bold';
            app.PlotLabel.Layout.Row = 1;
            app.PlotLabel.Layout.Column = 7;
            app.PlotLabel.Text = 'Plot:';

            % Create LineSlider
            app.LineSlider = uislider(app.GridLayout3);
            app.LineSlider.Limits = [1 50];
            app.LineSlider.ValueChangedFcn = createCallbackFcn(app, @LineSliderValueChanged, true);
            app.LineSlider.BusyAction = 'cancel';
            app.LineSlider.Layout.Row = 4;
            app.LineSlider.Layout.Column = [7 12];
            app.LineSlider.Value = 46.7804059722045;

            % Create PointEditFieldLabel
            app.PointEditFieldLabel = uilabel(app.GridLayout3);
            app.PointEditFieldLabel.HorizontalAlignment = 'right';
            app.PointEditFieldLabel.Layout.Row = 5;
            app.PointEditFieldLabel.Layout.Column = 7;
            app.PointEditFieldLabel.Text = 'Point';

            % Create PointEditField
            app.PointEditField = uieditfield(app.GridLayout3, 'numeric');
            app.PointEditField.ValueChangedFcn = createCallbackFcn(app, @PointEditFieldValueChanged, true);
            app.PointEditField.Layout.Row = 5;
            app.PointEditField.Layout.Column = 8;

            % Create imagesReconstrctionTab
            app.imagesReconstrctionTab = uitab(app.TabGroup);
            app.imagesReconstrctionTab.Title = 'Images Reconstrction';

            % Create UIAxesImg
            app.UIAxesImg = uiaxes(app.imagesReconstrctionTab);
            title(app.UIAxesImg, 'FLIm Images')
            app.UIAxesImg.Position = [663 263 530 392];

            % Create EditFieldLabel
            app.EditFieldLabel = uilabel(app.imagesReconstrctionTab);
            app.EditFieldLabel.HorizontalAlignment = 'right';
            app.EditFieldLabel.Position = [18 659 55 22];
            app.EditFieldLabel.Text = 'Edit Field';

            % Create EditField
            app.EditField = uieditfield(app.imagesReconstrctionTab, 'text');
            app.EditField.Position = [88 659 505 22];

            % Create DataFileSelectedLabel
            app.DataFileSelectedLabel = uilabel(app.UIFigure);
            app.DataFileSelectedLabel.HorizontalAlignment = 'right';
            app.DataFileSelectedLabel.FontWeight = 'bold';
            app.DataFileSelectedLabel.FontColor = [1 0 0];
            app.DataFileSelectedLabel.Position = [691 703 112 22];
            app.DataFileSelectedLabel.Text = 'Data File Selected:';

            % Create DataFileSelectedDisp
            app.DataFileSelectedDisp = uieditfield(app.UIFigure, 'text');
            app.DataFileSelectedDisp.Editable = 'off';
            app.DataFileSelectedDisp.Position = [806 702 438 22];

            % Create ContextMenu
            app.ContextMenu = uicontextmenu(app.UIFigure);

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = mainTriplexGui_exported

            runningApp = getRunningApp(app);

            % Check for running singleton app
            if isempty(runningApp)

                % Create UIFigure and components
                createComponents(app)

                % Register the app with App Designer
                registerApp(app, app.UIFigure)

                % Execute the startup function
                runStartupFcn(app, @startupFcn)
            else

                % Focus the running singleton app
                figure(runningApp.UIFigure)

                app = runningApp;
            end

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end