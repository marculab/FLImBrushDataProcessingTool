function out=alpha_up(M,K,init,eps)
% find the upper bound of alpha value
% input: init - initial guess of alpha value
%        eps - the upper bound for the condition number

if ~exist('init','var')||isempty(init)
    init=0.9;
end

if ~exist('eps','var')||isempty(eps)
    eps=1.0001;
end

L=@(alpha)Laguerre(M,K,alpha);
cond_val=@(alpha)cond(L(alpha)'*L(alpha))-eps;

out=fzero(cond_val,init);

end