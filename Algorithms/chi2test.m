function [h,pvalue,stat,cvalue]=chi2test(res,param)

T=size(res,1);
alpha=0.05;

if length(param)/2>=1
    
    for i_param=1:length(param)/2
       
        eval([genvarname(param{i_param*2-1}),'=param{i_param*2};']);
        
    end
        
end

if ~exist('vars','var')||isempty(vars)
%     vars = var(res(round(T*0.1):end,:)); % calculate variance for each decay
%     avg = mean(res(round(T*0.1):end,:)); % calculate mean residue for each decay
    vars=var(reshape(res(round(T*0.5):end,:),[],1));% Estimate variance from last 50% of data points
    
end

if ~exist('dof','var')||isempty(dof)
    
    dof=T;
    
end
% figure;plot(res)
% res = bsxfun(@minus,res,avg); % remove mean from residue
% figure;plot(res)
stat=sum(res.^2)./vars/dof;

pvalue=1-chi2cdf(stat*dof,dof);

cvalue=chi2inv(1-alpha,dof)/dof;

h=(alpha>=pvalue);

end