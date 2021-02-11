% code to check APD system waveform shift

%% Channel 1
rawData1 = Ch1DataObj.rawData;
figure
plot(rawData1(:,1:8)); % channel 1 shift is 1

%% Channel 2 
rawData2 = Ch2DataObj.rawData; 
figure
plot(rawData2(:,1:8)); % channel 2 shift is 1

%% Channel 3 
rawData3 = Ch3DataObj.rawData; 
%%
idx = 3738;
figure
plot(rawData3(:,4*(idx-1)+1:4*(idx+2))); % channel 3 shift is 1
CtrlV = Ch3DataObj.CtrlV;

%%
figure
plot(CtrlV)




