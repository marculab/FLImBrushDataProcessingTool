% BiExp test script
clear all
close all

bg = loadTRIPLEXBG('');
irf = loadTRIPLEXiIRF('');

[namelist, folder] = CreatInputFilelist('');
[name_cell, folder_cell] = creatJobQueue(namelist, folder,1);
num_of_batches = length(name_cell);
mkdir(folder_cell{1},'DeCon_BiExp');

for i=1:num_of_batches
    i
data = RawDataStrcutClass(name_cell(i),folder_cell(i),bg,irf);
timeStamp = get(data,'timeStamp');
[pathstr,name,ext] = fileparts(name_cell{i});
savepath = fullfile(folder_cell{i},'DeCon_BiExp',strcat(name_cell{i},'_BiExp_DeCon.mat'));

low = 0.05;
high = 0.6;
[out,non_decon_idx,full_data_size]= detectSaturation(data,low,high);

fullIdx = 1:full_data_size(1);
fullIdx=fullIdx';
decon_idx = setdiff(fullIdx,non_decon_idx);
% reduced_data = out;
% reduced_data(idx,:)=[];
% currentdata  = get(data,'rawdata');
% dif = currentdata-out;

%initialize final data matrix
LT1 = zeros(full_data_size(1),1);
LT2 = zeros(full_data_size(1),1);
Weight1 = zeros(full_data_size(1),2);
Weight2 = zeros(full_data_size(1),2);
INT1 = zeros(full_data_size(1),1);
INT2 = zeros(full_data_size(1),1);
Tau1 = zeros(full_data_size(1),2);
Tau2 = zeros(full_data_size(1),2);

if ~isempty(decon_idx)
    
curser_low = 170;
curser_high = 370;

[BG_removed, BG_scale_factor] = removeBG(data,curser_low,curser_high);


truncation_length=680;
peak_position =[443 1194 0 0];

datastruct = truncateData(data,truncation_length,peak_position);

dt = 0.08;

% data = datasruct{1}';

% mask = max(datasruct{1},[],2);
% mask=mask>0.05;
% data = data(:,mask);
% [A, T, avglife, intensity, fitt, raw] = multiexp_fit(data,dt,laser,2);
addpath(genpath('C:\Users\XZ\OneDrive\UC Davis\Projects\Programming\Matlab Code\trfs_triplex_gui'));

rawdatatruct1 = channeldata(datastruct{1}',irf{1},dt,1.5);
expstruct1=ExpModel(rawdatatruct1,2);
expstruct1.iIRF_align();
expstruct1.fit_exp();



rawdatatruct2 = channeldata(datastruct{2}',irf{2},dt,1.5);
expstruct2=ExpModel(rawdatatruct2,2);
expstruct2.iIRF_align();
expstruct2.fit_exp();

LT1(decon_idx,:) = expstruct1.LTs';
Weight1(decon_idx,:) = expstruct1.Weights';
INT1(decon_idx,:) = expstruct1.INTs';
LT2(decon_idx,:) = expstruct2.LTs';
Weight2(decon_idx,:) = expstruct2.Weights';
INT2(decon_idx,:) = expstruct2.INTs';
Tau1(decon_idx,:) = expstruct1.Taus';
Tau2(decon_idx,:) = expstruct2.Taus';
    
end
lifet_avg{1}=LT1;
lifet_avg{2}=LT2;
spec_int{1}=INT1;
spec_int{2}=INT2;

save(savepath,'lifet_avg','Weight1','Weight2','spec_int','Tau1','Tau2','timeStamp');

end