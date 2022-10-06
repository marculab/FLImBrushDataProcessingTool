clear
close all
clc

%% load file locations
load('N:\V2 Sacramento Database\Brain Necrosis\Data_Aggregate\matfile_adrs.mat')
idx = (T.Instrument=='FLImBrush');
fileNames = T.mat_file_adr;
fileNames = fileNames(idx);
%%
for i=9:317
    i
    temp = fileNames{i};
    [filepath,name,ext] = fileparts(temp);
    newName = [name '_12.5GS' ext];
    load(temp, 'EOP_H1G','EOP_H1S','SP_G','SP_S','calibrationObj')
    newFile = fullfile(filepath,newName);
    save(newFile,'EOP_H1G','EOP_H1S','SP_G','SP_S','calibrationObj','-append')
    
end