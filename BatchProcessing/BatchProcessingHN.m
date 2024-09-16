clear
close all
clc
addpath(genpath('..\Algorithms'))

%% Set Deconvolution Variables
root = 'D:\Data_100_Patient_Study';

APD1Path = '..\APDDetectorFile\UCD_MC\M00562768_laser_reflection.mat';
APD2Path = '..\APDDetectorFile\UCD_MC\M00549707_DCS.mat';
APD3Path = '..\APDDetectorFile\UCD_MC\M00549708_DASPI.mat';
APD4Path = '..\APDDetectorFile\UCD_MC\M00928868 DASPI.mat';
ChWidth = 680*0.08;
alpha1 = 0.916;
alpha2 = 0.916;
alpha3 = 0.916;
alpha4 = 0.965;
LGorder = 12;
 
%% Set up the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 26);

% Specify range and delimiter
opts.DataLines = [4, Inf];
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["Patient", "Run", "ScanContext", "Anatomy", "MCReferenceFrame", "DataChannelsUsed", "IncludeInAnalysis", "RawFLImDataFiletdms", "BackgroundFile", "DeconvolutionFile", "TextFile", "VideoPath", "NonMCCorrectedPositions", "RefinedPositions", "AnnotationFile", "ReferenceFrameWLI", "VarName17", "HistopathologyOrder", "Deconvolution", "ABPositionCorrection", "InstrumentSegmentation", "MotionCorrection", "HistopathologyCoregistration", "AnnotationMask", "ReprocessingDeconvolution", "OtherNotes"];
opts.VariableTypes = ["double", "double", "categorical", "categorical", "double", "categorical", "categorical", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "categorical", "categorical", "string", "categorical", "categorical", "string"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, ["RawFLImDataFiletdms", "BackgroundFile", "DeconvolutionFile", "TextFile", "VideoPath", "NonMCCorrectedPositions", "RefinedPositions", "AnnotationFile", "ReferenceFrameWLI", "VarName17", "HistopathologyOrder", "Deconvolution", "ABPositionCorrection", "HistopathologyCoregistration", "OtherNotes"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["ScanContext", "Anatomy", "DataChannelsUsed", "IncludeInAnalysis", "RawFLImDataFiletdms", "BackgroundFile", "DeconvolutionFile", "TextFile", "VideoPath", "NonMCCorrectedPositions", "RefinedPositions", "AnnotationFile", "ReferenceFrameWLI", "VarName17", "HistopathologyOrder", "Deconvolution", "ABPositionCorrection", "InstrumentSegmentation", "MotionCorrection", "HistopathologyCoregistration", "AnnotationMask", "ReprocessingDeconvolution", "OtherNotes"], "EmptyFieldRule", "auto");

% Import the data
HNAggregateInputDataRunLevel = readtable("C:\Users\MarcuAdmin\Documents\GitHub\FLImBrushDataProcessingTool\BatchProcessing\H&N Aggregate Input Data - Run_Level.csv", opts);

% Clear temporary variables
clear opts

%% Select Patients to be processed
RunLevelData = HNAggregateInputDataRunLevel;
T = RunLevelData(RunLevelData.Patient>204,:);


%% 
%for i = 1:height(T)
for i = 19:20
    ID = T.Patient(i);
    Run = T.Run(i);

    DataPath = T.RawFLImDataFiletdms(i);

    BGPath = T.BackgroundFile(i);

    disp(['Processing ' DataPath]);

    processRun_HN(root, DataPath, BGPath, APD1Path, APD2Path, APD3Path, APD4Path, alpha1, alpha2, alpha3, alpha4, ChWidth, LGorder);


end