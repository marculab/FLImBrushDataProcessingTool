clear
close all
clc

%% Import data from text file
% Script for importing data from the following text file:
%
%    filename: E:\MyGitRepo\FLImBrushDataProcessingTool\AdcancedFunctions\H&N Aggregate Input Data - Run_Level 20231106.csv
%
% Auto-generated by MATLAB on 06-Nov-2023 18:33:28

%% Set up the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 26);

% Specify range and delimiter
opts.DataLines = [2, Inf];
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["Patient", "Run", "ScanContext", "Anatomy", "MCReferenceFrame", "DataChannelsUsed", "IncludeInAnalysis", "RawFLImDataFiletdms", "BackgroundFile", "DeconvolutionFile", "TextFile", "VideoPath", "NonMCCorrectedPositions", "RefinedPositions", "AnnotationFile", "ReferenceFrameWLI", "VarName17", "HistopathologyOrder", "Deconvolution", "ABPositionCorrection", "InstrumentSegmentation", "MotionCorrection", "HistopathologyCoregistration", "AnnotationMask", "ReprocessingDeconvolution", "OtherNotes"];
opts.VariableTypes = ["double", "double", "categorical", "categorical", "double", "double", "categorical", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "categorical", "categorical", "string", "categorical", "categorical", "string"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, ["RawFLImDataFiletdms", "BackgroundFile", "DeconvolutionFile", "TextFile", "VideoPath", "NonMCCorrectedPositions", "RefinedPositions", "AnnotationFile", "ReferenceFrameWLI", "VarName17", "HistopathologyOrder", "Deconvolution", "ABPositionCorrection", "HistopathologyCoregistration", "OtherNotes"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["ScanContext", "Anatomy", "IncludeInAnalysis", "RawFLImDataFiletdms", "BackgroundFile", "DeconvolutionFile", "TextFile", "VideoPath", "NonMCCorrectedPositions", "RefinedPositions", "AnnotationFile", "ReferenceFrameWLI", "VarName17", "HistopathologyOrder", "Deconvolution", "ABPositionCorrection", "InstrumentSegmentation", "MotionCorrection", "HistopathologyCoregistration", "AnnotationMask", "ReprocessingDeconvolution", "OtherNotes"], "EmptyFieldRule", "auto");

% Import the data
HN = readtable("E:\MyGitRepo\FLImBrushDataProcessingTool\AdcancedFunctions\H&N Aggregate Input Data - Run_Level 20231106.csv", opts);
T = HN((HN.Patient>100),:);

%% Clear temporary variables
clear opts

%%
root = 'Y:\V2 Sacramento Database\Da Vinci Robot Study (100 patients)\Data_100_Patient_Study';
for i=313:size(T,1)
    i
    temp = T.DeconvolutionFile{i};
    if ~isempty(temp)
        matName = fullfile(root,temp)
        load(matName)
        %     temp(1) = 'W';
        % SNR1 = Ch1DataObj.SNR;
        %--------------------------channel 1-------------------------------
        Ch1DataObj.noise = std(Ch1DataObj.preProcessedData(end-200:end,:),1); % use last 200 points since data is upsampled
        WFMax = max(Ch1DataObj.preProcessedData);
        Ch1DataObj.SNR = 20*log10(WFMax./Ch1DataObj.noise)'; % covert to column vector
        % SNR11 = Ch1DataObj.SNR;
        % scatter(SNR1, SNR11)
        %--------------------------channel 2-------------------------------
        Ch2DataObj.noise = std(Ch2DataObj.preProcessedData(end-200:end,:),1); % use last 200 points since data is upsampled
        WFMax = max(Ch2DataObj.preProcessedData);
        Ch2DataObj.SNR = 20*log10(WFMax./Ch2DataObj.noise)'; % covert to column vector
        %--------------------------channel 3-------------------------------
        Ch3DataObj.noise = std(Ch3DataObj.preProcessedData(end-200:end,:),1); % use last 200 points since data is upsampled
        WFMax = max(Ch3DataObj.preProcessedData);
        Ch3DataObj.SNR = 20*log10(WFMax./Ch3DataObj.noise)'; % covert to column vector

        if ~exist('EOP_H1G','var')
            CH1WFGainCorr = get(Ch1DataObj,'GainCorrectedWF');
            CH2WFGainCorr = get(Ch2DataObj,'GainCorrectedWF');
            CH3WFGainCorr = get(Ch3DataObj,'GainCorrectedWF');
            %-------------------------------------------compute extended phasor-----------------------------------------------------
            WFExtended = [CH1WFGainCorr;CH2WFGainCorr;CH3WFGainCorr];
            numOfWF = size(WFExtended,2);
            EOP = zeros(numOfWF,1);
            for j = 1:numOfWF
                EOP(i) = ComputePhasor(WFExtended(:,i),0,1,0);
            end
            EOP_H1G = real(EOP);
            EOP_H1S = imag(EOP);
            %-------------------------------------------compute spectral phasor--------------------------------------------
            Spectrum = zeros(numOfWF,3);
            Spectrum(:,1) = sum(CH1WFGainCorr)';
            Spectrum(:,2) = sum(CH2WFGainCorr)';
            Spectrum(:,3) = sum(CH3WFGainCorr)';
            SP = zeros(numOfWF,1);
            for j = 1: numOfWF
                SP(i) = (ComputePhasor(Spectrum(i,:),0,1,0)+(1+1i))/2;
            end
            SP_G = real(SP);
            SP_S = imag(SP);
        end
        [filepath,name,ext] = fileparts(matName);
        newFileName = fullfile(filepath,[name '_new' ext]);
        save(newFileName, 'dataInfoObj','Ch1DataObj','Ch2DataObj','Ch3DataObj','EOP_H1G','EOP_H1S','SP_G','SP_S','-v7.3');
        delete(matName)
        status = movefile(newFileName,matName);
        saveDeconLite(Ch1DataObj,Ch2DataObj,Ch3DataObj,matName)
    end
end

