% test script
clear all
close all

bg = loadTRIPLEXBG('');
irf = loadTRIPLEXiIRF('');



data = RawDataStrcutClass('',bg,irf);

low = 0.05;
high = 0.6;
[out, idx] = detectSaturation(data,low,high);

% currentdata  = get(data,'rawdata');
% dif = currentdata-out;

curser_low = 170;
curser_high = 370;

BG_removed = removeBG(data,curser_low,curser_high);


truncation_length=680;
peak_position =[443 1194 1937 1937];
datasruct = truncateData(data,truncation_length,peak_position);

path('C:\Users\XZ\Desktop\BitBucket\trfs_triplex_gui\Algorithms\MultiExp Generic',path);
path('C:\Users\XZ\Desktop\BitBucket\trfs_triplex_gui\Algorithms',path);
dt = 0.08;
laser = irf{1};


data = datasruct{1}';

mask = max(datasruct{1},[],2);
mask=mask>0.05;
data = data(:,mask);
% [A, T, avglife, intensity, fitt, raw] = multiexp_fit(data,dt,laser,2);

rawdatatruct = channeldata(data,laser,dt,1.5);
expstruct=ExpModel(rawdatatruct,2);
expstruct.iIRF_align();
expstruct.fit_exp();