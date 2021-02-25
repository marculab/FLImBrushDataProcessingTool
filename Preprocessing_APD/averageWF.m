function [data_out,outlierFlag] = averageWF(data_in, numOfAvg)
numOfWF = size(data_in,2);
if mod(numOfWF,numOfAvg)
    error('Wrong number of average. Total number of waveform can not be divided by number of average');
else
    outlierFlag = zeros(numOfWF,1); % outlier flag
    numOfWFPostAvg = size(data_in,2)/numOfAvg;
    data_out = zeros(size(data_in,1),numOfWFPostAvg);
    
    for i = 1:numOfWFPostAvg
        temp = data_in(:,numOfAvg*(i-1)+1:numOfAvg*i);  
%         temp = alignWaveform_CFDNew(temp,3,0.4); % align waveforms
        WFMax = max(temp);
        outlierFlagTemp = isoutlier(WFMax);
        outlierFlag(numOfAvg*(i-1)+1:numOfAvg*i) = outlierFlagTemp;
%         sum(outlierFlagTemp)
        if sum(outlierFlagTemp)>1 % if more than 1 outlier, remove all WF
            outlierFlagTemp = ones(size(outlierFlagTemp));
            
        end
        temp(:,outlierFlagTemp==1) = [];
%         temp = alignWaveform_CFDNew(temp,3,0.4); % align waveforms
        avgTemp = nanmean(temp,2);
        data_out(:,i) = avgTemp;
    end
end
