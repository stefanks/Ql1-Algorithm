function step = ql1_istastep(g,tau,x,alphabar)
step= max(x-alphabar*(g+tau),0) - max(-x-alphabar*(-g+tau),0)-x;            
end

