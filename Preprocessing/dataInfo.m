classdef dataInfo < handle
    %UNTITLED4 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Access = public)
        dataFileNames;
        dataFolderName;
        apd1Name;
        apd1Folder;
        apd2Name;
        apd2Folder;
        apd3Name;
        apd3Folder;
        bgFileName;
        bgFileFolder;
        caliFileName;
        caliFileFolder;
        laserRepRate = NaN;
    end
    
    
    methods
        %constructor
        function obj = dataInfo(dataFileNamesIn, dataFolderNamesIn, ...
                apd1NameIn, apd1FolderIn, apd2NameIn, apd2FolderIn, ...
                apd3NameIn, apd3FolderIn, bgFileNameIn, bgFolderNameIn,...
                caliFileNameIn, caliFileFolderIn, laserRepRateIn)
            %UNTITLED4 Construct an instance of this class
            %   Detailed explanation goes here
            obj.dataFileNames = dataFileNamesIn;
            obj.dataFolderName = dataFolderNamesIn;
            obj.apd1Name = apd1NameIn;
            obj.apd1Folder = apd1FolderIn;
            obj.apd2Name = apd2NameIn;
            obj.apd2Folder = apd2FolderIn;
            obj.apd3Name = apd3NameIn;
            obj.apd3Folder = apd3FolderIn;
            obj.bgFileName = bgFileNameIn;
            obj.bgFileFolder = bgFolderNameIn;
            obj.caliFileName = caliFileNameIn;
            obj.caliFileFolder = caliFileFolderIn;
            obj.laserRepRate = laserRepRateIn;
        end
        
%         function pathToOutput = writeDataInfoToDisk(filename,pathToFolder)
%             %METHOD1 Summary of this method goes here
%             %   Detailed explanation goes here
%             outputArg = obj.Property1 + inputArg;
%         end
        
    end
end

