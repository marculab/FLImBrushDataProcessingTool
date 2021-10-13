function [namelist, folder] = creatInputFileList(rootpath)

namelist = [];
folder = [];
while isempty(namelist)
    
    option = questdlg('Do you want to load a folder or specific files?',...
        'What to Load?', 'Folder', 'Files','Cancel', 'Folder');
    % folder option
    if strcmp(option,'Folder')
        
        folderpath = uigetdir(rootpath, 'Please Select the Raw Data Folder from TRIPLEX');
        if ~folderpath
            f=warndlg('No folder is selected, Please select again.');
            %             nameListing = [];
            waitfor(f);
            return
        else
            % 1. Read the folder contents
            dirListing = dir(folderpath);
            nameListing = {dirListing.name};
            % 1. Find and exclude folders within the folder
            find_Dirs = [dirListing.isdir];
            nameListing(find_Dirs) = [];
            %             name = nameListing;
            %             folder = folderpath;
        end
        
    end
    % files option
    if strcmp(option,'Files')
        [filename, folderpath] = uigetfile('*.tdms*', 'Please Select the Raw Data from TRIPLEX', 'Multiselect', 'on',rootpath);
        if ~folderpath
            f=warndlg('No data file is selected, Please select again.');
            waitfor(f);
            return
        end
        
        if ~iscell(filename)
            filename = {filename};
        end
        
        nameListing = filename;
        % 1. Find and exclude folders within the folder
%         find_Dirs = isdir(nameListing);
%         nameListing(find_Dirs) = [];
        %         name = nameListing;
        %         folder = folderpath;
    end
    
    if strcmp(option,'Cancel')
        nameListing = [];
    end
    
    
    
%     if ~isempty(nameListing)
%                 
%         % 2. Find and exclude filenames with extensions
%         find_dots = strfind(nameListing, '.');
%         nameListing(~cellfun('isempty', find_dots)) = [];
%         
%         % 3. Find and exclude lifetimes files
%         find_dots = strfind(nameListing, 'lifetimes');
%         nameListing(~cellfun('isempty', find_dots)) = [];
%         
%         % 4. Find and exclude header files
%         find_dots = strfind(nameListing, 'header');
%         nameListing(~cellfun('isempty', find_dots)) = [];
%         
%     end
    
    if isempty(nameListing)
        
        f=warndlg('No data has been found, Please select again.');
        waitfor(f);
        return
    else
        namelist = nameListing';
%         num_of_files = length(namelist);
        folder = cell(size(namelist));
        folder(:) = {folderpath};
    end
    
end
return
