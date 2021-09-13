function L=Laguerre(M,K,alpha)
% M is length
% K is the order
%  alpha is alpha
L = zeros(M,K);
for  j=1:K
    if j==1
        % Zero-order Laguerre function
        for i=1:M
            L(i,j)=sqrt(alpha^(i-1)*(1-alpha)); % discrete laguerre is non-negative, thus i-1
        end %for i
    else
        %  All other Laguerre functions
        for i=1:M
            if i==1
                L(i,j)=sqrt(alpha)*L(i,j-1); %changed from k to k-1(from the paper)
            else
                L(i,j)=sqrt(alpha)* (L(i-1,j)+L(i,j-1))-L(i-1,j-1);
            end
        end
    end
end