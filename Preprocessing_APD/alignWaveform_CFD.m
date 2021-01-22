function [data_out] = alignWaveform_CFD(data_in,f,peakMaxLocs)
% in-input matrix with waveform as columns
% data with peak aligned 
switch nargin
    case 1
        f = 0.5;
        peakMaxLocs = size(data_in,1);
    case 2
        peakMaxLocs = size(data_in,1);
    case 3
        if peakMaxLocs > size(data_in,1)
            warning('Peak max location exceed data size, set peakMaxLocs to data size.')
            peakMaxLocs = size(data_in,1);
        end
    otherwise
        error('Wrong number of Imputs')
end

x = 1:size(data_in,1);
% why did I do interplation?
%xx = 0.2:0.2:size(data_in,1);
%data_in = interp1(x,data_in,xx);
data_in(isnan(data_in)) = eps;
data_out = data_in; % make copy
[maxV,maxIdxs]=max(data_in(1:peakMaxLocs,:)); % find max index
maxIdxs = max(maxIdxs);
risingEdge = data_in(1:maxIdxs,:);
risingEdge(1:maxIdxs-40,:)=0;
%--------------------------------------------------------------------------
% figure
% plot(risingEdge(:,1:100:end))
% title('Rising Edge')
%--------------------------------------------------------------------------
risingEdge = bsxfun(@minus,risingEdge,f*maxV);
risingEdge(risingEdge>0)=inf;
%--------------------------------------------------------------------------
% figure
% plot(risingEdge(:,1:100:end))
% title('Rising Edge - half amplitude')
%--------------------------------------------------------------------------
risingEdge = abs(risingEdge);
% figure
% plot(risingEdge(:,1:100:end))
% title('Abs of Rising Edge - half amplitude')

%--------------------------------------------------------------------------

[~,CFDIdx] = min(risingEdge);
% histogram(CFDIdx)
shift = mode(CFDIdx)-CFDIdx; % calculate shift
for i = 1:length(shift) %loop through all columns
    data_out(:,i) = circshift(data_in(:,i),shift(i));
end
%data_out = downsample(data_out,5);