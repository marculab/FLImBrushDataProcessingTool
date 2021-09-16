function y = triexp_model(x,a1,a2,a3,t1,t2,t3,L)
%global X
decay = a1.*exp(-x./t1)+a2.*exp(-x./t2)+a3.*exp(-x./t3);
%decay = decay./max(decay);
y = filter(L,1,decay);
%y = conv(X,a1.*exp(-x./t1)+(1-a1).*exp(-x./t2));
% y = y(1:length(x));
% y = y./max(y);
end