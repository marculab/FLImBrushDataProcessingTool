function IntuitiveSurgicalDataAggregate = importIntuitiveDataAggregate(filename, dataLines)
%IMPORTFILE Import data from a text file
%  INTUITIVESURGICALDATAAGGREGATE = IMPORTFILE(FILENAME) reads data from
%  text file FILENAME for the default selection.  Returns the data as a
%  table.
%
%  INTUITIVESURGICALDATAAGGREGATE = IMPORTFILE(FILE, DATALINES) reads
%  data for the specified row interval(s) of text file FILENAME. Specify
%  DATALINES as a positive scalar integer or a N-by-2 array of positive
%  scalar integers for dis-contiguous row intervals.
%
%  Example:
%  IntuitiveSurgicalDataAggregate = importfile("E:\MyGitRepo\FLImBrushDataProcessingTool\BatchProcessing\Intuitive Surgical Data Aggregate.csv", [2, Inf]);
%
%  See also READTABLE.
%
% Auto-generated by MATLAB on 07-Apr-2023 16:32:00

%% Input handling

% If dataLines is not specified, define defaults
if nargin < 2
    dataLines = [2, Inf];
end

%% Set up the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 37);

% Specify range and delimiter
opts.DataLines = dataLines;
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["DataGroup", "Run", "SpecimenSiteID", "SpecimenOrSiteDescription", "TimeScanStart", "TimeScanEnd", "MidpointOfScanTime", "ScanContext", "MotionCorrection", "IncludeInAnalysis", "FiberID", "DaVinciUsed", "SealingType", "BurstingLeftRightmmHg", "OtherNotes", "EnergyApplied", "EnergyInstrument", "PowerDelivered", "Duration", "TimeEnergyWasApplied", "DeltaBetweenCauteryTimeAndScanTime", "DataChannelsUsed", "BGToUse", "ReferenceFrame", "HistopathologyRegistered", "AnnotationTool", "DeconvolutionPerformed", "PositionCorrected", "VarName29", "RawDataFile", "BackgroundFile", "DeconvolutionFile", "TextFile", "VideoPath", "NonMCCorrectedPositions", "RefinedPositions", "AnnotationFile"];
opts.VariableTypes = ["categorical", "double", "string", "categorical", "string", "string", "string", "categorical", "categorical", "categorical", "double", "categorical", "categorical", "string", "string", "categorical", "double", "categorical", "double", "string", "string", "categorical", "double", "double", "categorical", "categorical", "categorical", "categorical", "string", "string", "categorical", "categorical", "categorical", "string", "categorical", "categorical", "categorical"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, ["SpecimenSiteID", "TimeScanStart", "TimeScanEnd", "MidpointOfScanTime", "BurstingLeftRightmmHg", "OtherNotes", "TimeEnergyWasApplied", "DeltaBetweenCauteryTimeAndScanTime", "VarName29", "RawDataFile", "VideoPath"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["DataGroup", "SpecimenSiteID", "SpecimenOrSiteDescription", "TimeScanStart", "TimeScanEnd", "MidpointOfScanTime", "ScanContext", "MotionCorrection", "IncludeInAnalysis", "DaVinciUsed", "SealingType", "BurstingLeftRightmmHg", "OtherNotes", "EnergyApplied", "PowerDelivered", "TimeEnergyWasApplied", "DeltaBetweenCauteryTimeAndScanTime", "DataChannelsUsed", "HistopathologyRegistered", "AnnotationTool", "DeconvolutionPerformed", "PositionCorrected", "VarName29", "RawDataFile", "BackgroundFile", "DeconvolutionFile", "TextFile", "VideoPath", "NonMCCorrectedPositions", "RefinedPositions", "AnnotationFile"], "EmptyFieldRule", "auto");
opts = setvaropts(opts, ["FiberID", "EnergyInstrument", "Duration", "BGToUse"], "TrimNonNumeric", true);
opts = setvaropts(opts, ["FiberID", "EnergyInstrument", "Duration", "BGToUse"], "ThousandsSeparator", ",");

% Import the data
IntuitiveSurgicalDataAggregate = readtable(filename, opts);

end