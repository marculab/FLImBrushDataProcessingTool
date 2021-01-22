function [h,pvalue,stat,cvalue]=lbqtest(res,param)
% param is the cell of parameters: lags, alpha, dof

T=size(res,1);

lags=min(20,T-1);
alpha=0.05;
dof=lags;

if length(param)/2>=1
    
    for i_param=1:length(param)/2
       
        eval([genvarname(param{i_param*2-1}),'=param{i_param*2};']);
        
    end
        
end

if lags>T-1
    
    lags=min(20,T-1);
    warndlg('Lags must not exceed the length of the data minus one; default value is used.');
    
end

if dof>lags
    
    dof=lags;
    warndlg('Degrees of freedom must not exceed corresponding lag; default value is used.');
    
end

% Autocorrelation
nFFT=2^(nextpow2(size(res,1))+1);
F=fft(bsxfun(@minus,res,mean(res,1)),nFFT);
F=F.*conj(F);
acf=ifft(F);
acf=acf(1:(lags+1),:);%Retain non-negative lags
acf=bsxfun(@rdivide,acf,acf(1,:));%Normalize
acf=real(acf(2:end,:));%Strip off ACF at lag 0

% Compute Q-statistics to the largest lag; keep only those requested:
idx=(T-(1:lags))';
stat=T*(T+2)*cumsum(bsxfun(@rdivide,(acf.^2),idx));
stat=stat(lags,:);

% Compute p-values:
pvalue=1-chi2cdf(stat,dof);

% Compute critical values
cvalue=chi2inv(1-alpha,dof);

% Perform the test:
h=(alpha>=pvalue);

end