function ret = ql1_fValue(g,b,tau,x)
ax =  g + b;
ret = 1/2 * x' * ax - b' * x + sum(tau.*abs(x));
end
