function HNAggregateInputData20240822 = importHNRunData(workbookFile, sheetName, dataLines)
%IMPORTFILE Import data from a spreadsheet
%  HNAGGREGATEINPUTDATA20240822 = IMPORTFILE(FILE) reads data from the
%  first worksheet in the Microsoft Excel spreadsheet file named FILE.
%  Returns the data as a table.
%
%  HNAGGREGATEINPUTDATA20240822 = IMPORTFILE(FILE, SHEET) reads from the
%  specified worksheet.
%
%  HNAGGREGATEINPUTDATA20240822 = IMPORTFILE(FILE, SHEET, DATALINES)
%  reads from the specified worksheet for the specified row interval(s).
%  Specify DATALINES as a positive scalar integer or a N-by-2 array of
%  positive scalar integers for dis-contiguous row intervals.
%
%  Example:
%  HNAggregateInputData20240822 = importfile("E:\MyGitRepo\FLImBrushDataProcessingTool\AdcancedFunctions\H&N Aggregate Input Data 20240822.xlsx", "Run_Level", [1, Inf]);
%
%  See also READTABLE.
%
% Auto-generated by MATLAB on 22-Aug-2024 21:33:49

%% Input handling

% If no sheet is specified, read first sheet
if nargin == 1 || isempty(sheetName)
    sheetName = 1;
end

% If row start and end points are not specified, define defaults
if nargin <= 2
    dataLines = [1, Inf];
end

%% Set up the Import Options and import the data
opts = spreadsheetImportOptions("NumVariables", 26);

% Specify sheet and range
opts.Sheet = sheetName;
opts.DataRange = dataLines(1, :);

% Specify column names and types
opts.VariableNames = ["Patient", "Run", "ScanContext", "Anatomy", "MCReferenceFrame", "DataChannelsUsed", "IncludeInAnalysis", "RawFLImDataFiletdms", "BackgroundFile", "DeconvolutionFile", "TextFile", "VideoPath", "NonMCCorrectedPositions", "RefinedPositions", "AnnotationFile", "ReferenceFrameWLI", "VarName17", "HistopathologyOrder", "Deconvolution", "ABPositionCorrection", "InstrumentSegmentation", "MotionCorrection", "HistopathologyCoregistration", "AnnotationMask", "ReprocessingDeconvolution", "OtherNotes"];
opts.VariableTypes = ["double", "double", "categorical", "categorical", "double", "categorical", "categorical", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "categorical", "categorical", "string", "categorical", "categorical", "string"];

% Specify variable properties
opts = setvaropts(opts, ["RawFLImDataFiletdms", "BackgroundFile", "DeconvolutionFile", "TextFile", "VideoPath", "NonMCCorrectedPositions", "RefinedPositions", "AnnotationFile", "ReferenceFrameWLI", "VarName17", "HistopathologyOrder", "Deconvolution", "ABPositionCorrection", "HistopathologyCoregistration", "OtherNotes"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["ScanContext", "Anatomy", "DataChannelsUsed", "IncludeInAnalysis", "RawFLImDataFiletdms", "BackgroundFile", "DeconvolutionFile", "TextFile", "VideoPath", "NonMCCorrectedPositions", "RefinedPositions", "AnnotationFile", "ReferenceFrameWLI", "VarName17", "HistopathologyOrder", "Deconvolution", "ABPositionCorrection", "InstrumentSegmentation", "MotionCorrection", "HistopathologyCoregistration", "AnnotationMask", "ReprocessingDeconvolution", "OtherNotes"], "EmptyFieldRule", "auto");

% Import the data
HNAggregateInputData20240822 = readtable(workbookFile, opts, "UseExcel", false);

for idx = 2:size(dataLines, 1)
    opts.DataRange = dataLines(idx, :);
    tb = readtable(workbookFile, opts, "UseExcel", false);
    HNAggregateInputData20240822 = [HNAggregateInputData20240822; tb]; %#ok<AGROW>
end

end