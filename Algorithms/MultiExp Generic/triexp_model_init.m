function y = triexp_model_init(x,a1,a2,a3,t1,t2,t3)
%global X
y = a1.*exp(-x./t1)+a2.*exp(-x./t2)+a3.*exp(-x./t3);
% y = y./max(y);
%y = filter(X,1,decay);
%y = conv(X,a1.*exp(-x./t1)+(1-a1).*exp(-x./t2));
%y = y(1:length(x));
%y = y./max(y);
end 