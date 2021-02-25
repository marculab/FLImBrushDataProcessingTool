function out=findTruncationL(K,alpha,init,eps)
% find the upper bound of alpha value
% input: init - initial guess of alpha value
%        eps - the upper bound for the condition number

if ~exist('init','var')||isempty(init)
    init=272;
end

if ~exist('eps','var')||isempty(eps)
    eps=0.0001;
end

trunL = round(init*0.50):round(init*10);
trunL = trunL';
cond_val = zeros(size(trunL));
for i = 1:length(trunL)
L = Laguerre(trunL(i),K,alpha);
cond_val(i) = cond(L'*L)-1;
end

idx = find(cond_val<eps, 1 );

out=trunL(idx);

end