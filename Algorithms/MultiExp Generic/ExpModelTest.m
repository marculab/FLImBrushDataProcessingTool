clear
close all
clc

%% 
load('..\..\TestData\ExpTestData.mat')
%% 
N = 1000;
tau1 = rand(N,1)*4;
tau2 = rand(N,1)*4+4;
t = 0:1:length(laser)-1;
t = t*0.1;
decay = zeros(length(laser),N);
for i = 1:N
decay(:,i) = 0.5*exp(-t/tau1(i))+0.5*exp(-t/tau2(i));
end
figure;plot(decay)
WF = filter(laser,1,decay);
figure;plot(WF)

O = ExpModel(2, WF,laser, 0.1, '', '', '',[]);
%%
runDecon(O)
%%
idx = ceil(rand*N)
fit = get(O,'fit',idx);
figure
plot(WF(:,idx))
hold on
plot(fit)
hold off
%%
figure
scatter(O.LTs_formula,O.LTs_formula-O.LTs_decay);
