clear
close all
clc

%% 
load('..\..\TestData\ExpTestData.mat')
%% 
O = ExpModel(4, specTempCh1,laser, 0.2, '', '', '');
%%
runDecon(O)

fit = get(O,'fit');

%%
idx = 20
figure
plot(specTempCh1(:,idx))
hold on
plot(fit(:,idx))
hold off
