function ret = ql1_phi(g,tau,x)
ret = g + tau.*sign(x);
ret(x==0)=0;
end