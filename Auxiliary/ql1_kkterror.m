function finalerror = ql1_kkterror(g,tau,x)

finalerror = 0;

for i = 1: size(x,1)
    u=max(x(i),0);
    if g(i)+tau(i)<0
        uerror = -g(i)-tau(i);
    else
        uerror = (g(i)+tau(i))*(u/(u+1));
    end
    
    finalerror = max(finalerror, uerror);
    
    v=max(-x(i),0);
    if -g(i)+tau(i)<0
        verror = g(i)-tau(i);
    else
        verror = (-g(i)+tau(i))*(v/(v+1));
    end
    
    finalerror = max(finalerror, verror);
end
% 
% finalerror2 = 0;
% 
% for i = 1: size(x,1)
%     u=max(x(i),0);
%     v=max(-x(i),0);
%     if g(i)>tau(i)
%         finalerror2 = max(finalerror2, (g(i)+tau(i))*(u/(u+1)));
%         finalerror2 = max(finalerror2, g(i)-tau(i));
%     elseif g(i)<-tau(i)
%         finalerror2 = max(finalerror2, -g(i)-tau(i));
%         finalerror2 = max(finalerror2,(-g(i)+tau(i))*(v/(v+1)));  
%     else            
%         finalerror2 = max(finalerror2, (g(i)+tau(i))*(u/(u+1)));
%         finalerror2 = max(finalerror2,(-g(i)+tau(i))*(v/(v+1))); 
%     end
% end
% 
% if norm(finalerror2-finalerror)>1e-16
%     fprintf('ERORRRRRRRRRR');
%     finalerror= 'error';
% end

end

