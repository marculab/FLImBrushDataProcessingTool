function [data_out,dataVariationFlag] = averageWF(data_in, numOfAvg)
numOfWF = size(data_in,2);
if mod(numOfWF,numOfAvg)
    error('Wrong number of average. Total number of waveform can not be divided by number of average');
else
    numOfWFPostAvg = size(data_in,2)/numOfAvg;
    data_out = zeros(size(data_in,1),numOfWFPostAvg);
    dataVariationFlag = zeros(numOfWFPostAvg,1);
    for i = 1:numOfWFPostAvg
        temp = data_in(:,numOfAvg*(i-1)+1:numOfAvg*i);
        avgTemp = nanmean(temp,2);
        [MaxV,MaxIdx] = max(avgTemp);
        data_out(:,i) = avgTemp;
        data_std = nanstd(temp,[],2);
        dataMaxvar = data_std(MaxIdx)/MaxV;
        if dataMaxvar>0.01
            dataVariationFlag(i) = dataMaxvar;
        end
    end
end