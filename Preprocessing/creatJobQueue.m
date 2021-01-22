function [name_cell, folder_cell] = creatJobQueue (namelist, folderlist,limit)
% this is the function that divide all the files into pipline job list, in
% put is root path and limit for a single decovolution. out put is a cell
% array of matix. The 1st column is the file file name, 2nd column is the
% file folder
% namelist: column vector of all file names
% folderlist: column vector of all folder path

% remove duplicated folder path
% [unique_folder, ia, ic]= unique(folderlist);
% find total number of data folder
% num_of_folder = size(unique_folder);
% limit = 1*e9;

%colume cell array to store size of each data file


full_name = fullfile(folderlist,namelist);
temp_dir = dir(full_name{1});
singel_file_size_info=temp_dir.bytes;
if singel_file_size_info>limit
    name_cell = namelist;
    folder_cell = folderlist;
else
    
size_info = cell(size(namelist));

num_of_file = size(namelist,1);

for i=1:num_of_file
    temp_dir = dir(full_name{i});
    size_info{i}=temp_dir.bytes;
end
cum_size = cumsum(cell2mat(size_info));
% size_total = sum(cell2mat(size_info));
clear temp_dir

%initilize index and counter
I = 1;
II=1;
counter = 1;


while II<=num_of_file
        
    temp = (cum_size - limit);
    
    I = max(find(temp<0))
    
    name_cell{counter} = namelist(II:I);
    folder_cell{counter} = folderlist(II:I);
    cum_size = cum_size-cum_size(I);
    II=I+1;
    
    counter = counter+1;
    
end
end



