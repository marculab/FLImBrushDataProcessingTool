function y = quadriexp_model_init(x,a1,a2,a3,a4,t1,t2,t3,t4)

y = a1.*exp(-x./t1)+a2.*exp(-x./t2)+a3.*exp(-x./t3)+a4.*exp(-x./t4);

end