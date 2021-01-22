function [h,autoCorrCurve] = autoCorrTest(res, threshold)

T=size(res,2);
decayLength = size(res,1);
h = zeros(1,T);
autoCorrCurve = cell(T,1);
for i = 1:T
    autoCorrCurveTemp = xcorr(res(:,i), 'coeff'); %calculate cross correlation, normalized
    autoCorrCurve{i} = autoCorrCurveTemp(decayLength:end); %take the second half of xcorr
    res_temp = abs(res(:,i)); % absolute value of residue
    res_temp = res_temp>0.05; % point idx that larger than 0.05
    resPercent = double(sum(res_temp))/decayLength; % bad residue presentage
    
    temp = abs(autoCorrCurveTemp(decayLength:end));
    temp = temp>threshold;
    autoCorrPercent = double(sum(temp))/decayLength; % bad autocorelation presentage
    if (autoCorrPercent<0.5)&&(resPercent<0.1) % percentage of the data point to check the goodness of fit
        h(i) = 1;
    end
end