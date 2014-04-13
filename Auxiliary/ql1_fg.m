function [ f ,g ] = ql1_fg( Ax, b,x )
ax =  Ax(x);
f = 1/2 * x' * ax - b' * x ;
g=ax-b;
end