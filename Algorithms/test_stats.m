function stats = test_stats(spec,spec_fitted,dt,bw,varargin)
% test_stats(raw waveform, fitted wavefrom, time resolution, bandwidth)
% function to evaluate goodness of fit of deconvolution
% (Laguerre, multi-exponential or any other fitting algorithm)
% to 

threshold = 2/sqrt(size(spec,1));% threshod for auto-correlaton
if nargin > 4
    param = varargin{1};
else
    param = {[],[],threshold}; %2/sqrt(dacayLength) is the default value for the bound of max xcorr.
end

res = spec - spec_fitted; %residue

lbq = struct('h',[],'pvalue',[],'stat',[],'cvalue',[]);
chi2 = struct('h',[],'pvalue',[],'stat',[],'cvalue',[]);
autoCorr = struct('h',[],'autoCorrCurve',[]);
stats = struct('lbq',lbq,'chi2',chi2,'autoCorr',autoCorr);

[stats.autoCorr.h,stats.autoCorr.autoCorrCurve] = autoCorrTest(res,param{3});



% Downsampling residue
n_t=round(1/bw/dt);
% n_t=1;
[~,n_1]=max(spec_fitted);
n_1=quantile(n_1,0.95)+size(res,1)*0.05;
if n_1>size(res,1)-n_t*20
    n_1=1;
end

res=res(n_1:n_t:end,:);

% [stats.lbq.h,stats.lbq.pvalue,stats.lbq.stat,stats.lbq.cvalue]=lbqtest(res,param{1});
% [stats.chi2.h,stats.chi2.pvalue,stats.chi2.stat,stats.chi2.cvalue]=chi2test(res,param{2});
end