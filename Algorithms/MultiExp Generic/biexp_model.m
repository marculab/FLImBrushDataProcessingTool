function y = biexp_model(x,a1,t1,t2,L)
%global X
decay = a1.*exp(-x./t1)+(1-a1).*exp(-x./t2);
%decay = decay./max(decay);
y = filter(L,1,decay);
%y = conv(X,a1.*exp(-x./t1)+(1-a1).*exp(-x./t2));
y = y(1:length(x));
y = y./max(y);
end 