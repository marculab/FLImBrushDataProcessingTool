% code to determin FLImBRUSH decon parameter
%% 
clear all
close all
clc

%% load in iRF
load('../APDDetectorFile/M00549708_DASPI.mat')
%%
dt = irfdt;
LtPpix = 4; % ns
t = (1:5000)*dt;
decay =exp(-t/LtPpix);
figure
plot(t,decay)
xlabel('ns')
ylabel('Amplituede')
title('PpIX decay , tau = 16 ns')
%% calculate digitizer resolution
ENOB = 8.0; %50 Ohm impedience, NI 5162
Vres = 2/power(2,ENOB);
%% convolve sloest irf with decay
iIRF = irfUpSampled(:,210);
figure('Position',[200 200 600 600*0.75])
plot(iIRF)
rawData = filter(iIRF,1,decay);
rawData = rawData/max(rawData);
[~,MaxIdx] = max(rawData);
IdxTemp= find(rawData(MaxIdx:end)<Vres, 1 );
figure('Position',[200 200 600 600*0.75])
plot(t,rawData)
hold on
yline(Vres,'-.r','Digitizer resolution ENOB = 8','LineWidth',1)
xline(t(MaxIdx),'-.b',num2str(t(MaxIdx)))
xline(t(IdxTemp+MaxIdx-1),'-.b',num2str(t(IdxTemp+MaxIdx-1)))
hold off
text(65,0.7,['Tail Window(80%) = ',num2str(IdxTemp*dt),' ns']);
text(65,0.6,['Full Window = ',num2str((IdxTemp*dt)*5/4),' ns']);
xlabel('Time (ns)')
ylabel('Amplituede')
title('4 ns decay convolved with slowest irf , tau = 4 ns')
 %% convolve fastest irf with decay
iIRF = irfUpSampled(:,end);
figure('Position',[200 200 600 600*0.75])
plot(iIRF)
rawData = filter(iIRF,1,decay);
rawData = rawData/max(rawData);
[~,MaxIdx] = max(rawData);
IdxTemp= find(rawData(MaxIdx:end)<Vres, 1 );
figure('Position',[200 200 600 600*0.75])
plot(t,rawData)
hold on
yline(Vres,'-.r','Digitizer resolution ENOB = 8','LineWidth',1)
xline(t(MaxIdx),'-.b',num2str(t(MaxIdx)))
xline(t(IdxTemp+MaxIdx-1),'-.b',num2str(t(IdxTemp+MaxIdx-1)))
hold off
text(65,0.7,['Tail Window(80%) = ',num2str(IdxTemp*dt),' ns']);
text(65,0.6,['Full Window = ',num2str((IdxTemp*dt)*5/4),' ns']);
xlabel('Time (ns)')
ylabel('Amplituede')
title('4 ns decay convolved with fastest irf , tau = 4 ns')

%% determine alpha value
FW = ceil((IdxTemp*dt)*5/4);
L = ceil(FW/0.4)
out1=alpha_up(L,12,0.9,1.01);
out2=alpha_up(L,12,0.9);
disp([out1 out2])
