clear
close all
clc

%% load file locations
load('N:\V2 Sacramento Database\Brain Necrosis\Data_Aggregate\matfile_adrs.mat')
%%
idx = (T.Instrument=='FLImBrush');
fileNames = T.mat_file_adr;
fileNames = fileNames(idx);

for i=1:2
    
    temp = fileNames{i};
%     temp(1) = 'W';
    ReprocessFB(temp)
    
end

