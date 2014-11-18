function ret = ql1_phitilde(g,tau,x,alphabar)
ret= (x-max(abs(x-alphabar*g)-alphabar*tau,0).*sign(x-alphabar*g))/alphabar;
ret(x==0)=0;
end