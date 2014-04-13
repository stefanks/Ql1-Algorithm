function ret = ql1_vTrue(AvectorProduct,b,tau,x)
ret=AvectorProduct(x)-b+tau.*sign(x);
ret((abs(ret)<=tau & x==0))=0;
ret(x==0)=ret(x==0)-sign(ret(x==0)).*tau(x==0);
end