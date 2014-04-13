function fullHistory=alg_sub_RecordMV(fullHistory,i,g,b,tau,x,optimalityMeasure, steptype)
    fullHistory.fValues(i)= ql1_fValue(g,b,tau,x);
    fullHistory.sparsity(i) = size(find(x==0),1);
    
    % type of step. Numeric.
    fullHistory.MVTypes(i)=steptype;
    
    % This is the provided optimality measure
    fullHistory.optimalityMeasures(i) = optimalityMeasure(g,b,tau,x);
    fullHistory.normV(i) = norm(ql1_v(g,tau,x),Inf);
    fullHistory.kkterror(i) = ql1_kkterror(g,tau,x);
end
