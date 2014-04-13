function ret = ql1_fValueTrue(Ax,b,tau,x)
ret = 1/2 * x' * Ax(x) - b' * x + sum(tau.*abs(x));
end