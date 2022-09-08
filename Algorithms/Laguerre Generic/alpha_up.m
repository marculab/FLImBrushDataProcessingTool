function out=alpha_up(M,K,init,eps)
% M is data length
% K is Laguerre Order
% find the upper bound of alpha value
% input: init - initial guess of alpha value
%        eps - the upper bound for the condition number

if ~exist('init','var')||isempty(init)
    init=0.9;
end

if ~exist('eps','var')||isempty(eps)
    eps=1.01; % updated by Xiangnan 20210907
%     eps=1.000000001; % updated by Xiangnan 20220829 for testing
end

L=@(alpha)Laguerre(M,K,alpha);
cond_val=@(alpha)cond(L(alpha)'*L(alpha))-eps;

out=fzero(cond_val,init);

end