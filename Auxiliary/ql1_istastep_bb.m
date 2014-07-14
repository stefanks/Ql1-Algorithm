function [xF, gF, numMV,fullHistory,prevfValuesForIstaBB,xPrevOutput] = ql1_istastep_bb(problem,g,x,alphabar,s,As,optimalityMeasure,accuracy, numMV ,maxMV, fullHistory,nargoutGlobal,prevfValuesForIstaBB, MforIstaBB, xi,xPrevOutput,outputLevel,stallEpsilon)
%% Constants
nu = 2;          % suggested >1. Used 2

%% Compute intial stepsize
bbstepsize=min(s'*s/(s'*As),1/stallEpsilon);

%% Actual line search procedure



if isfield(problem, 'B')
    
    while 1
        xF=max(x-bbstepsize*(g+problem.tau),0) - max(-x-bbstepsize*(-g+problem.tau),0);
        Bmy = problem.B(xF)-problem.y;
        if 1/2 *(Bmy'*Bmy) -(1/2)*(problem.y'*problem.y) + (1/2)*problem.gamma*(xF'*xF) + sum(problem.tau.*abs(xF)) <= max(prevfValuesForIstaBB) - xi/(bbstepsize*2) *norm(xF-x)^2;
            break;
        end
        if bbstepsize<=alphabar
            break;
        end
        bbstepsize=max(bbstepsize/nu,alphabar);
    end
    gF = problem.Bt(Bmy)+problem.gamma*xF;
    if size(prevfValuesForIstaBB,1)==MforIstaBB+1
        prevfValuesForIstaBB(1:MforIstaBB)=prevfValuesForIstaBB(2:MforIstaBB+1);
        prevfValuesForIstaBB(MforIstaBB+1) = ql1_fValue(gF,problem.b,problem.tau,xF);
    else
        prevfValuesForIstaBB = [prevfValuesForIstaBB;ql1_fValue(gF,problem.b,problem.tau,xF)];
    end
else
    
    while 1
        xF=max(x-bbstepsize*(g+problem.tau),0) - max(-x-bbstepsize*(-g+problem.tau),0);
        gF=problem.Ax(xF)-problem.b;
        numMV = numMV+1;
        if (nargoutGlobal>=2),fullHistory=alg_sub_RecordMV(fullHistory, numMV,gF,problem.b,problem.tau,xF,optimalityMeasure,1);end
        if (outputLevel >=3), xPrevOutput=alg_sub_OutputX(gF,problem.b,problem.tau,xF,optimalityMeasure,numMV,'*xBB',(nargoutGlobal >= 3), xPrevOutput); end
        if numMV>=maxMV
            break;
        end
        if ql1_fValue(gF,problem.b,problem.tau,xF) <= max(prevfValuesForIstaBB) - xi/(bbstepsize*2) *norm(xF-x)^2;
            break;
        end
        if optimalityMeasure(gF,problem.b,problem.tau,xF) <= accuracy
            break;
        end
        if bbstepsize<=alphabar
            break;
        end
        bbstepsize=max(bbstepsize/nu,alphabar);
    end
    if size(prevfValuesForIstaBB,1)==MforIstaBB+1
        prevfValuesForIstaBB(1:MforIstaBB)=prevfValuesForIstaBB(2:MforIstaBB+1);
        prevfValuesForIstaBB(MforIstaBB+1) = ql1_fValue(gF,problem.b,problem.tau,xF);
    else
        prevfValuesForIstaBB = [prevfValuesForIstaBB;ql1_fValue(gF,problem.b,problem.tau,xF)];
    end
    
end


end