function [data_out,varargout] = alignWaveform_CFDNew(data_in, t_rising, dt, f_in, range_in)
% in-input matrix with waveform as columns
% f is fraction
% t_rising is the rising time of the signal
% dt is sampling interval
% range is the range of data used for alignment
switch nargin
    case 3
        f = 0.5;
        range = 1:size(data_in,1);
    case 4
        f = f_in;
        range = 1:size(data_in,1);
    case 5
        f = f_in;
        range = range_in;
    otherwise
        error('Wrong number of Inputs')
end
% remove NaN
data_in_raw = data_in; % make copy of raw data
data_in(isnan(data_in)) = eps;

% remove signal DC value
DC_range = round(size(data_in,1)*0.1);
DC = mean(mean(data_in(1:DC_range,:)));
data_in = data_in-DC;
% figure;plot(data_in(:,1:50:end));
%---------------truncate data----------------------
data_in = data_in(range,:);
% delay positive part
t_d = ceil(t_rising/dt*0.5);
% data_in(end-t_d-1:end,:) = data_in(end-t_d-1:end,:)*0;
dataDelayed = circshift(data_in,t_d,1);
%
dataInvert = -data_in*f;
% figure
% plot(dataDelayed(:,1:100:end))
% hold on
% plot(dataInvert(:,1:100:end))

dataSum = dataDelayed+dataInvert;
% figure;plot(dataSum(:,1:100:end))
dataSum(1:t_d+1,:) = dataSum(1:t_d+1,:)*0;
% figure;plot(dataSum(:,4)); hold on; plot(dataInvert(:,4));plot(dataDelayed(:,4));hold off
[~,minIdx] = min(dataSum);
minIdx = floor(mode(minIdx))-1;
if minIdx<1
    minIdx = 1;
end

[~,maxIdx] = max(dataSum);
maxIdx = round(mean(maxIdx))+1;
maxIdx = min(maxIdx,size(dataSum,1));
% plotIdx = 100;
% figure;plot(dataDelayed(:,plotIdx));hold on;plot(dataSum(:,plotIdx));hold off; 
% pause(1);
% close(gcf);
risingEdge = dataSum(minIdx:maxIdx,:);
% figure;plot(risingEdge(:,1:500:end))
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
% figure;plot(abs(risingEdge))
[~,zeroCrossingIdxMin] = min(abs(risingEdge));
% parfor i = 1:size(risingEdge,2) %loop through all columns
%     data_out(:,i) = circshift(data_in_raw(:,i),shift(i));
% end
% zeroCrossingIdxMin = zeroCrossingIdxMin';
shift = mode(zeroCrossingIdxMin)-zeroCrossingIdxMin; % calculate shift
data_out = data_in_raw;
parfor i = 1:length(shift) %loop through all columns
    data_out(:,i) = circshift(data_in_raw(:,i),shift(i));
end
varargout{1} = shift;
% figure;plot(data_out(:,1:100:end))
%data_out = downsample(data_out,5); 