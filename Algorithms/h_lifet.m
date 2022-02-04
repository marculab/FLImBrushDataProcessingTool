function [LT,INT]=h_lifet(h,dt,type,p)
% dt - time resolution
% type - {'average','1/e','quantile'}
% p - (1) type=='1/e', p is the 1/e location
%   - (2) type=='quantile', p is the quantile

t_bin=0:size(h,1)-1;


if ~exist('type','var')||isempty(type)
    
    type='average';
    
end

switch type
    
    case 'average'
        
        lifet_sampl=(0.5+t_bin)*h./sum(h);%mean
        lifet_sampl = lifet_sampl';
        
    case '1/e'
        
        if ~exist('p','var')||isempty(p)
    
            p=1/exp(1);
    
        end

        
        h_surf=bsxfun(@rdivide,h,max(h));
        t1=sum(h_surf>p)';
        
        t1(t1<=0)=1;
        t1(t1>=size(h,1))=t1(t1>=size(h,1))-1;
        
        t1=t1+(p-diag(h_surf(t1,:)))./(diag(h_surf(t1+1,:))-diag(h_surf(t1,:)));
        lifet_sampl=t1-0.5;
        
    case 'quantile'
        
        if ~exist('p','var')||isempty(p)
    
            p=1-1/exp(1);
    
        end
        
        h_cumsum=cumsum(bsxfun(@rdivide,h,sum(h)));
        t1=sum(h_cumsum<=p)';
        
        t1(t1<=0)=1;
        t1(t1>=size(h,1))=t1(t1>=size(h,1))-1;
        
        t1=t1+(p-diag(h_cumsum(t1,:)))./(diag(h_cumsum(t1+1,:))-diag(h_cumsum(t1,:)));
        lifet_sampl=t1-0.5;
        
        
end

LT=lifet_sampl*dt;
INT=sum(h);