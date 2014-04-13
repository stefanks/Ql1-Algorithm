function outX = alg_sub_OutputX(g,b,tau,x,optimalityMeasure,numMV,varname,recordAllHistory, varargin)
if nargin==8
    fprintf(' numMV xType    OptMeasure        fValue        num x<0            num x==0           num x>0\n');
    fprintf('%6d %5s %13.3e %13.3e %6i             %6i             %6i\n' ,numMV, varname, optimalityMeasure(g,b,tau,x), ql1_fValue(g,b,tau,x), size(find(x<0),1),size(find(x==0),1),size(find(x>0) ,1));
    outX =x;
else
    xHistory=varargin{1};
    xPrevOutput=xHistory(:,end);
    if norm(x-xPrevOutput)<1e-16 && size(find(sign(x)-sign(xPrevOutput)==0),1) == size(x,1)
        outX=xHistory;
    else
        numnegadd = size(find(xPrevOutput>=0 & x<0),1);
        numzeroadd = size(find(xPrevOutput~=0 & x==0),1);
        numposadd = size(find(xPrevOutput<=0 & x>0),1);
        numnegremoved = size(find(xPrevOutput<0 & x>=0),1);
        numzeroremoved = size(find(xPrevOutput==0 & x~=0),1);
        numposremoved = size(find(xPrevOutput>0 & x<=0),1);
        fprintf('%6d %5s %13.3e %13.3e %6i' ,numMV, varname, optimalityMeasure(g,b,tau,x), ql1_fValue(g,b,tau,x), size(find(x<0),1));
        if numnegadd>0 || numnegremoved>0
            fprintf('(');
            if numnegadd>0
                fprintf('+%4i',numnegadd);
            else
                fprintf('     ');
            end
            if numnegremoved>0
                fprintf('-%4i',numnegremoved);
            else
                fprintf('     ');
            end
            fprintf(')');
        else
            fprintf('            ');
        end
        fprintf(' %6i', size(find(x==0),1));
        
        if numzeroadd>0 || numzeroremoved>0
            fprintf('(');
            if numzeroadd>0
                fprintf('+%4i',numzeroadd);
            else
                fprintf('     ');
            end
            if numzeroremoved>0
                fprintf('-%4i',numzeroremoved);
            else
                fprintf('     ');
            end
            fprintf(')');
        else
            fprintf('            ');
        end
        
        fprintf(' %6i', size(find(x>0),1));
        
        if numposadd>0 || numposremoved>0
            fprintf('(');
            if numposadd>0
                fprintf('+%4i',numposadd);
            else
                fprintf('     ');
            end
            if numposremoved>0
                fprintf('-%4i',numposremoved);
            else
                fprintf('     ');
            end
            fprintf(')');
        else
            fprintf('            ');
        end
        
        fprintf('\n');
               
        if recordAllHistory && nargin==9
            outX =[xHistory x];
        else
            outX =x;
        end
    end
end
end
