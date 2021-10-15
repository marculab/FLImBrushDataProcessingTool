function data_out = averageGianOrV(data_in, numOfAvg)
    numOfDataPostAvg = length(data_in)/numOfAvg;
    data_out = zeros(numOfDataPostAvg,1); % column vector
    for i = 1:numOfDataPostAvg
        temp = data_in((i-1)*numOfAvg+1:i*numOfAvg);
        tempMean = mean(temp,'omitnan');
        if all(temp==tempMean) % check whether all data within averaging window are the same, of so average data , if not throw error
            data_out(i) = tempMean;
        else
            error('Differnt Gain or Control V within averaging window, possible corruption of data!');
        end
    end
    
