function ret = ql1_omega(g,tau,x) 
ret = g-tau.*sign(g);
ret(x~=0)=0;
ret(abs(g)-tau<=0)=0;
end