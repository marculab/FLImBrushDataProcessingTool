function LCs = nnls3c(spec,basis,iIRF)

vv = filter(iIRF,1,basis);
D_mat=conv2(eye(size(vv,1)),[1,-3,3,-1],'valid')'*basis;
D=D_mat;
H=vv'*vv;
H_chol=chol(inv(H));
C=H_chol*D';
l1=H_chol*vv';

lam=zeros(size(D,1),size(spec,2));
parfor i=1:size(spec,2)
    d=l1*spec(:,i);
    lam(:,i)=lsqnonneg(C,d);
end
LCs=(vv'*vv)\(vv'*spec-D'*lam);
LCs = LCs./(ones(size(LCs,1),1)*sum(vv*LCs));
end