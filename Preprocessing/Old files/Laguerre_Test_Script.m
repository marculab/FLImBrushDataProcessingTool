% BiExp test script
clear all
close all

BG = loadTRIPLEXBG('');
[irf,~] = loadTRIPLEXiIRF('');

[namelist, folder] = CreatInputFilelist('');
[name_cell, folder_cell] = creatJobQueue(namelist, folder,1);
num_of_batches = length(name_cell);
mkdir(folder_cell{1},'DeCon_Laguerre');

%% loop for each line
for i=1:num_of_batches
    i
    lifet_avg = cell(4,1);
    spec_int = cell(4,1);
    alpha_out = cell(4,1);
    SNR = cell(4,1);
    test = cell(4,1);
    Laguerre_coeffs = cell(4,1);
    
data = RawDataStrcutClass(name_cell(i),folder_cell(i),BG,irf);

timeStamp = get(data,'timeStamp');
[pathstr,name,ext] = fileparts(name_cell{i});
savepath = fullfile(folder_cell{i},'DeCon_Laguerre',strcat(name_cell{i},'_Laguerre_DeCon.mat'));
% get gain information
[gain_factor, gain_list] = gainCorrection(data);

low = 0.07;
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
LT3 = zeros(full_data_size(1),1);
LT4 = zeros(full_data_size(1),1);
% Weight1 = zeros(full_data_size(1),2);
% Weight2 = zeros(full_data_size(1),2);
INT1 = zeros(full_data_size(1),1);
INT2 = zeros(full_data_size(1),1);
INT3 = zeros(full_data_size(1),1);
INT4 = zeros(full_data_size(1),1);

SNR1 = ones(full_data_size(1),1)*50;
SNR2 = ones(full_data_size(1),1)*50;
SNR3 = ones(full_data_size(1),1)*50;
SNR4 = ones(full_data_size(1),1)*50;

BG_scale_factor = zeros(full_data_size(1),1);;
% Tau1 = zeros(full_data_size(1),2);
% Tau2 = zeros(full_data_size(1),2);
CH1_LCs = zeros(full_data_size(1),12);
CH2_LCs = zeros(full_data_size(1),12);

if ~isempty(decon_idx)
    
curser_low = 1000;
curser_high = 1200;

[BG_removed, BG_scale_factor] = removeBG(data,curser_low,curser_high);


truncation_length=680;
peak_position =[515 1281 0 0];

datastruct = truncateData(data,truncation_length,peak_position);

dt = 0.08;

% data = datasruct{1}';

% mask = max(datasruct{1},[],2);
% mask=mask>0.05;
% data = data(:,mask);
% [A, T, avglife, intensity, fitt, raw] = multiexp_fit(data,dt,laser,2);
addpath(genpath('C:\Users\XZ\OneDrive\UC Davis\Projects\Programming\Matlab Code\trfs_triplex_gui'));

rawdatatruct1 = channeldata(datastruct{1}',irf{1},dt,1.5);
LaguerreModel_CH1=LaguerreModel(rawdatatruct1);
LaguerreModel_CH1.iIRF_align();
LaguerreModel_CH1.estimate_laguerre();

% h=figure;
% plot(rawdatatruct1.data);
% hold on
% 
% plot(get(LaguerreModel_CH1,'fit'));

rawdatatruct2 = channeldata(datastruct{2}',irf{2},dt,1.5);
LaguerreModel_CH2=LaguerreModel(rawdatatruct2);
LaguerreModel_CH2.iIRF_align();
LaguerreModel_CH2.estimate_laguerre();

LT1(decon_idx,:) = LaguerreModel_CH1.LTs';
% Weight1(decon_idx,:) = expsLaguerreModel_CH1.Weights';
INT1(decon_idx,:) = LaguerreModel_CH1.INTs';
LT2(decon_idx,:) = LaguerreModel_CH2.LTs';
% Weight2(decon_idx,:) = expstruct2.Weights';
INT2(decon_idx,:) = LaguerreModel_CH2.INTs';
% Tau1(decon_idx,:) = expsLaguerreModel_CH1.Taus';
% Tau2(decon_idx,:) = expstruct2.Taus';
CH1_LCs(decon_idx,:) = LaguerreModel_CH1.LCs';
CH2_LCs(decon_idx,:) = LaguerreModel_CH2.LCs';

alpha_out{1} = get(LaguerreModel_CH1,'alpha');
alpha_out{2} = get(LaguerreModel_CH2,'alpha');
alpha_out{3} = [];
alpha_out{4} = [];



end

lifet_avg{1}=LT1;
lifet_avg{2}=LT2;
lifet_avg{3}=LT3;
lifet_avg{4}=LT4;

spec_int{1}=INT1;
spec_int{2}=INT2;
spec_int{3}=INT3;
spec_int{4}=INT4;

SNR{1} = SNR1;
SNR{2} = SNR2;
SNR{3} = SNR3;
SNR{4} = SNR4;


Laguerre_coeffs{1} = CH1_LCs;
Laguerre_coeffs{2} = CH2_LCs;
% Laguerre_coeffs{3} = [];
% Laguerre_coeffs{4} = [];

% save to .mat file
save(savepath,'lifet_avg','spec_int','alpha_out','Laguerre_coeffs',...
    'timeStamp','SNR','test','decon_idx','full_data_size','BG',...
    'BG_scale_factor','gain_factor','gain_list');

end
disp('Finished');