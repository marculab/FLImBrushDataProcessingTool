function y = monoexp_model(x,a,t,L)
%global X
decay = a.*exp(-x./t);
%decay = decay./max(decay);
y = filter(L,1,decay);
%y = conv(X,a1.*exp(-x./t1)+(1-a1).*exp(-x./t2));
y = y(1:length(x));
y = y./max(y);
end 