function [xF, gF, numMV,fullHistory,prevfValuesForIstaBB,xPrevOutput] = ql1_istastep_bb(Ax, b,g,tau,x,alphabar,s,As,optimalityMeasure,accuracy, numMV ,maxMV, fullHistory,nargoutGlobal,prevfValuesForIstaBB, MforIstaBB, xi,xPrevOutput,outputLevel,stallEpsilon)
%% Constants
nu = 2;          % suggested >1. Used 2

%% Compute intial stepsize
bbstepsize=min(s'*s/(s'*As),1/stallEpsilon);
%% Actual line search procedure
while 1
    xF=max(x-bbstepsize*(g+tau),0) - max(-x-bbstepsize*(-g+tau),0);
    gF=Ax(xF)-b;
    numMV = numMV+1;
    if (nargoutGlobal>=2),fullHistory=alg_sub_RecordMV(fullHistory, numMV,gF,b,tau,xF,optimalityMeasure,1);end 
    if (outputLevel >=3), xPrevOutput=alg_sub_OutputX(gF,b,tau,xF,optimalityMeasure,numMV,'*xBB',(nargoutGlobal >= 3), xPrevOutput); end
    if numMV>=maxMV
        break;
    end
    if ql1_fValue(gF,b,tau,xF) <= max(prevfValuesForIstaBB) - xi/(bbstepsize*2) *norm(xF-x)^2;
        break;
    end
    if optimalityMeasure(gF,b,tau,xF) <= accuracy
        break;
    end
    if bbstepsize<=alphabar
        break;
    end
    bbstepsize=max(bbstepsize/nu,alphabar);
end
if size(prevfValuesForIstaBB,1)==MforIstaBB+1
    prevfValuesForIstaBB(1:MforIstaBB)=prevfValuesForIstaBB(2:MforIstaBB+1);
    prevfValuesForIstaBB(MforIstaBB+1) = ql1_fValue(gF,b,tau,xF);
else
    prevfValuesForIstaBB = [prevfValuesForIstaBB;ql1_fValue(gF,b,tau,xF)];
end
end