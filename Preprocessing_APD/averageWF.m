function [data_out,outlierFlag] = averageWF(data_in, numOfAvg)
% average waveform
% if there are only one outlier, it is removed, the rest of WF are averaged
% if there are more than one outlier, all waveform are removed
% output will be NaN if all data in the average window are NaN
numOfWF = size(data_in,2);
if mod(numOfWF,numOfAvg)
    error('Wrong number of average. Total number of waveform can not be divided by number of average');
else
    outlierFlag = zeros(numOfWF,1); % outlier flag
    numOfWFPostAvg = size(data_in,2)/numOfAvg;
    data_out = zeros(size(data_in,1),numOfWFPostAvg);
    
    for i = 1:numOfWFPostAvg
        temp = data_in(:,numOfAvg*(i-1)+1:numOfAvg*i);  
%         temp = alignWaveform_CFDNew(temp_in,2.8, 0.4);
%         figure
%         tiledlayout(2,1)
%         nexttile
%         plot(temp_in,'*-');
%         nexttile
%         plot(temp,'*-');
%         temp = alignWaveform_CFDNew(temp,3,0.4); % align waveforms
        WFMax = max(temp); % find waveform maximum
        outlierFlagTemp = isoutlier(WFMax); % find outlier, NaN will not be treated as outlier
        outlierFlag(numOfAvg*(i-1)+1:numOfAvg*i) = outlierFlagTemp;
%         sum(outlierFlagTemp)
        if sum(outlierFlagTemp)>1 % if more than 1 outlier, remove all WF
            outlierFlagTemp = ones(size(outlierFlagTemp));
        end
        temp(:,outlierFlagTemp==1) = [];
        avgTemp = nanmean(temp,2);
        data_out(:,i) = avgTemp;
    end
end
