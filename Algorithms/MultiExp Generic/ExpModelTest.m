clear
close all
clc

%% 
load('..\..\TestData\ExpTestData.mat')
%% 
specTempCh1(:,2) = zeros(size(specTempCh1(:,2)));
laser_m = repmat(laser,1, size(specTempCh1,2));
laser_m(:,2) = zeros(size(laser_m(:,2)));
O = ExpModel(3, specTempCh1,laser, 0.2, '', '', '',[]);
%%
runDecon(O)
%%
idx = 10

fit = get(O,'fit',idx);

%%
figure
plot(specTempCh1(:,idx))
hold on
plot(fit)
hold off

figure
scatter(O.LTs_formula,O.LTs_formula-O.LTs_decay);
