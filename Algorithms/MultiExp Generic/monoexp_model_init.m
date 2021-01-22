function y = monoexp_model_init(x,a,t)
%global X
y = a.*exp(-x./t);
% y = y./max(y);
%y = filter(X,1,decay);
%y = conv(X,a1.*exp(-x./t1)+(1-a1).*exp(-x./t2));
%y = y(1:length(x));
%y = y./max(y);
end 