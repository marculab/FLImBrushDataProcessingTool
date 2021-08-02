function [data_out] = alignWaveform_CFDNew(data_in,t_rising, dt,f)
% in-input matrix with waveform as columns
% f is fraction
% t_rising is the rising time of the signal
% dt is sampling interval
% data with peak aligned
switch nargin
    case 3
        f = 0.6;
    case 4
        
    otherwise
        error('Wrong number of Inputs')
end
% remove NaN
data_in_raw = data_in; % make copy of raw data
data_in(isnan(data_in)) = eps;

% remove signal DC value
DC = mean(data_in(1:75));
data_in = data_in-DC;
% figure;plot(data_in(:,1:50:end));

% delay positive part
t_d = ceil(t_rising/dt*0.5);
dataDelayed = circshift(data_in,t_d,1);
%
dataInvert = -data_in*f;

%
dataSum = dataDelayed+dataInvert;
% figure;plot(dataSum(:,4)); hold on; plot(dataInvert(:,4));plot(dataDelayed(:,4));hold off
[~,minIdx] = min(dataInvert);
minIdx = mode(minIdx)-1;
if minIdx<1
    minIdx = 1;
end

% 
dataSum = dataDelayed+dataInvert;
% figure;plot(dataSum(:,4)); hold on; plot(dataInvert(:,4));plot(dataDelayed(:,4));hold off
[~,minIdx] = min(dataInvert);
minIdx = round(mean(minIdx))-1;
if minIdx<1
    minIdx = 1;
end

[~,maxIdx] = max(dataDelayed);
maxIdx = round(mean(maxIdx))+1;

% plotIdx = 100;
% figure;plot(dataDelayed(:,plotIdx));hold on;plot(dataSum(:,plotIdx));hold off; 
% pause(1);
% close(gcf);
risingEdge = dataSum(minIdx:maxIdx,:);
% zeroCrossingIdx =zeros(size(risingEdge,2),1);
% 
% for i = 1:size(risingEdge,2)
%     risingEdgeShifted = circshift(risingEdge(:,i),-1);
%     idxTemp = find((risingEdgeShifted.*risingEdge(:,i))<=0,1);
%     if isempty(idxTemp)
%         idxTemp = 0;
%     end
%     zeroCrossingIdx(i) = idxTemp+1;
% end
% shift = mode(zeroCrossingIdx)-zeroCrossingIdx; % calculate shift
[~,zeroCrossingIdxMin] = min(abs(risingEdge));
% zeroCrossingIdxMin = zeroCrossingIdxMin';
shift = mode(zeroCrossingIdxMin)-zeroCrossingIdxMin; % calculate shift
data_out = data_in_raw;
for i = 1:length(shift) %loop through all columns
    data_out(:,i) = circshift(data_in_raw(:,i),shift(i));
end
%data_out = downsample(data_out,5); 