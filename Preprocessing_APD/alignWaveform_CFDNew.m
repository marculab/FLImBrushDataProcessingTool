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

[~,maxIdx] = max(dataDelayed);
maxIdx = mode(maxIdx)+1;

if maxIdx-minIdx <8
    maxIdx = minIdx+8;
end
% plotIdx = 100;
% figure;plot(dataDelayed(:,plotIdx));hold on;plot(dataSum(:,plotIdx));hold off;
% pause(1);
% close(gcf);
risingEdge = dataSum(minIdx:maxIdx,:);
% find index of 1st positive after zero crossing and its value
[zeroCrossingValue,zeroCrossingIdx] = min(abs(risingEdge));
zeroCrossingIdx(zeroCrossingIdx<2)=2; % avoid out of bound error
zeroCrossingIdx(zeroCrossingIdx>(size(risingEdge,1)-1))=size(risingEdge,1)-1; % avoid out of bound error

prezeroCrossingValue = zeros(size(risingEdge,2),1);
postzeroCrossingValue = zeros(size(risingEdge,2),1);
% zeroCrossingValue = zeros(size(risingEdge,2),1);
for i = 1:size(risingEdge,2)
    prezeroCrossingValue(i) = risingEdge(zeroCrossingIdx(i)-1,i);
    postzeroCrossingValue(i) = risingEdge(zeroCrossingIdx(i)+1,i);
end
% for i = 1:size(risingEdge,2)
%     risingEdgeShifted = circshift(risingEdge(:,i),-1);
%     idxTemp = find((risingEdgeShifted.*risingEdge(:,i))<=0,1);
%     if isempty(idxTemp)
%         idxTemp = 0;
%     end
%     zeroCrossingIdx(i) = idxTemp+1;
%     prezeroCrossingValue(i) = risingEdge(idxTemp,i);
%     zeroCrossingValue(i) = risingEdge(idxTemp+1,i);
%     postzeroCrossingValue(i) = risingEdge(idxTemp+2,i);
% end
outLierFlag= isoutlier(zeroCrossingValue,'gesd');

if sum(outLierFlag) % if there is an outlier
%     figure;plot(risingEdge,'*-')
    outLierIdx = find(outLierFlag==1); % find outlier index
    outLierValue = zeroCrossingValue(outLierIdx); % get outlier value
    meanTemp = mean(zeroCrossingValue(~outLierFlag)); % get the mean of the rest
    if outLierValue > meanTemp % if outlier value is higher
        newValue = prezeroCrossingValue(outLierIdx); % use value from previous data
        newCrossingValue = zeroCrossingValue;
        newCrossingValue(outLierIdx) = newValue; % check whether it is still outlier
        if sum(isoutlier(newCrossingValue,'gesd')) % still a outlier
            data_in_raw(:,outLierIdx) = NaN(size(data_in_raw(:,outLierIdx)));
        else % not an outlier, reduce shift by 1
            zeroCrossingIdx(outLierIdx) = zeroCrossingIdx(outLierIdx)-1;
        end
    else % if outlier is lower
        newValue = postzeroCrossingValue(outLierIdx); % use value from previous data
        newCrossingValue = zeroCrossingValue;
        newCrossingValue(outLierIdx) = newValue; % check whether it is still outlier
        if sum(isoutlier(newCrossingValue,'gesd')) % still a outlier
            data_in_raw(:,outLierIdx) = NaN(size(data_in_raw(:,outLierIdx)));
        else % not an outlier, increase shift by 1
            zeroCrossingIdx(outLierIdx) = zeroCrossingIdx(outLierIdx)+1;
        end
    end
    % shift = mode(zeroCrossingIdx)-zeroCrossingIdx; % calculate shift
    % find cloest index near zero crossing
end
    shift = mode(zeroCrossingIdx)-zeroCrossingIdx;
    data_out = data_in_raw;
    for i = 1:length(shift) %loop through all columns
        if ~sum(isnan(data_in_raw(:,i)))
            data_out(:,i) = circshift(data_in_raw(:,i),shift(i));
        end
    end

%data_out = downsample(data_out,5);