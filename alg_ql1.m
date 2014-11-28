function varargout = alg_ql1(problem,varargin)
%% Universal Algorithm
%
%   A method for finding an optimal solution to
%       min (1/2)*x'*A*x - b'*x + norm(tau.*x,1)
%
%   Must have these elements:
%       problem.Ax               - handle for matrix vector product
%       problem.b                - vector b
%       problem.tau              - vector or scalar tau
%
%   Optional:
%
%       problem.normA            - An upper bound on the norm of matrix A
%                                  (default = 1e6)
%       opts.optimalityMeasure   - Function handle for measuring optimality
%                                  (default = scaled ISTA step)
%                                  Other choices: norm(ql1_v), ql1_kkterror, ql1_fValue
%       opts.accuracy            - Accuracy that opt measure must reach
%                                  (default = -1, see guessOptimal)
%       opts.maxA                - Maximum number of Ax calls allowed
%                                  (default = 10000)
%       opts.c                   - Decrease parameter for CG step
%                                  (default = 1e-8)
%       opts.outputLevel         - How much output to display
%                                  0 - none
%                                  1 - Beginning and end
%                                  2 - Outer loop progress (default)
%                                  3 - Inner loop progress
%       opts.x_0                 - Starting point (default=all zeros)
%       opts.separate            - Separate into free vs zero steps
%                                  0 no
%                                  1 yes (default)
%       opts.guessOptimal        - Stop if guessing optimality
%                                  0 no
%                                  1 yes (default)
%
%
%     (NOT RECOMMENDED) If the problem is of the elastic net form
%       min (1/2)*norm(B*x-y)^2 + gamma*norm(x)^2 + norm(tau.*x,1)
%     Can provide these additional elements:
%       problem.B                - handle for matrix vector product B*x
%       problem.Bt               - handle for matrix vector product B'*x
%       problem.y                - vector y
%       problem.gamma            - vector gamma
%
%     Providing these values MAY speed up the line search, but only if
%     computation time of B(x) is faster than A(x). Also, optimality
%     is not checked in linesearch with this option.
%
%
%
%
%
%   Outputs:
%     1st argument:
%       x                        - Best solution found
%       g                        - Gradient of the smooth part of F at x
%       fValue                   - F evaluated at x
%       optMeasure               - Optimality measure at x
%       numZeros                 - Sparsity of x
%       numA                     - Number of Ax calls
%       numB                     - Number of Bx or B'x calls
%       numOuterIterations       - Number of outer iterations (CG cycles)
%       algStatus                - Final algorithm status. Can be:
%                                  'optimal' - found solution to required accuracy
%                                  'maxA'   - work limit reached
%                                  'stall' - no progress during a full outer iteration
%                                  'guessOptimal' - guess that opt solution found
%                                  'unbounded' - problem is unbounded
%     2nd argument: (optional and expensive! Only for debugging. Computed after each Ax call)
%       fValues                  - F value of each iterate
%       sparsity                 - Sparsity of each iterate
%       MVTypes                  - What were Ax calls used for:
%                                  0 - Initial g computation
%                                  1 - Used for first order step
%                                  2 - Used for relaxation step
%                                  3 - Used for CG step
%       optimalityMeasures       - Optimality measures
%       normV                    - Norms of minimum norm subgradient
%       kkterror                 - KKT errors of reformulated problem
%       CGmvcount                - Ax calls for each CG cycle
%       CGglobalMVlink           - Current MV count after each CG cycle
%       reasonForCGstop          - Why CG cycles were stopped:
%                                  'optimal' - Found optimal solution
%                                  'g balance' - Gradient balance triggered
%                                  'F no dec' - Not enough decrease in F
%                                  'q inc'     - Too much accuracy asked for
%                                  'stall d'   - Too much accuracy asked for
%                                  'stall dAd' - Too much accuracy asked for
%     3rd argument: (optional and expensive! Only for debugging. Computed after each Ax call)
%       xPrevOutput              - A sequence of all iterates

%% Add necessary auxilary files
warning('off','MATLAB:dispatcher:pathWarning');
addpath('Auxiliary');

%% Default Constants
stallingEpsilon=1e-24;     %CHECK THIS
M=5; % Line search memory parameter
xi=0.005; % Line search sufficient decrease parameter
nu=2; % Line search step decrease parameter

%% Read inputs
if nargin>1
    opts=varargin{1};
else
    opts=struct;
end
if isfield(problem, 'normA')
    alphabar = 1/problem.normA;
else
    alphabar = 1/(1e6);
end
if ~isfield( opts, 'optimalityMeasure' )
    optimalityMeasure = @(g,b,tau,x) norm(max(x/alphabar-(g+tau),0) - max(-x/alphabar-(-g+tau),0) - x/alphabar,Inf);
else
    optimalityMeasure  =  opts.optimalityMeasure;
end
if  ~isfield( opts, 'accuracy' )
    accuracy = -1;
else
    accuracy  =  opts.accuracy;
end
if  ~isfield( opts, 'maxA' )
    maxA = 20000;
else
    maxA  = opts.maxA;
end
if  ~isfield( opts, 'c' )
    c = 1e-8;
else
    c  = opts.c;
end
if  ~isfield( opts, 'outputLevel' )
    outputLevel = 2;
else
    outputLevel  = opts.outputLevel;
end
if ~isfield( opts, 'x_0' )
    x_0 = zeros(size(problem.b));
else
    x_0  =  opts.x_0;
end
if ~isfield( opts, 'separate' )
    separate = 1;
else
    separate  =  opts.separate;
end
if ~isfield( opts, 'guessOptimal' )
    guessOptimal = 1;
else
    guessOptimal  =  opts.guessOptimal;
end

%% Load problem parameters from input
Ax = problem.Ax;
b = problem.b;
tau = problem.tau;

%% Configure the gradient balance
gb=@(g,tau,x)norm(ql1_omega(g,tau,x))^2 <= norm(ql1_phitilde(g,tau,x,alphabar))^2;

%% Initialize variables for output
numA = 0;
numB=0;
numOuterIterations = 0;
fullHistory=struct;

if nargout >= 2
    % history for each CG iteration
    fullHistory.reasonForCGstop = cell(maxA+2,1);
    fullHistory.CGmvcount = zeros(maxA+2,1);
    fullHistory.CGglobalMVlink=zeros(maxA+2,1);
    
    % history for each MV product
    fullHistory.fValues = zeros(maxA+2,1);
    fullHistory.MVTypes = zeros(maxA+2,1);
    fullHistory.optimalityMeasures = zeros(maxA+2,1);
    fullHistory.sparsity = zeros(maxA+2,1);
    fullHistory.normV= zeros(maxA+2,1);
    fullHistory.kkterror= zeros(maxA+2,1);
end

%% Compute needed starting information
x=x_0;
g=Ax(x)-b;
numA = numA+1;
if (nargout >=2),fullHistory=alg_sub_RecordMV(fullHistory, numA,g,b,tau,x,optimalityMeasure,0);end

%% Output: initial
if outputLevel>=1
    fprintf('-----Starting alg_ql1-----  ');
    fprintf('\n');
    xPrevOutput = alg_sub_OutputX(g,b,tau,x,optimalityMeasure,numA,numB,'x',(nargout >= 3));
else
    xPrevOutput = [];
end

%% Needed for ISTA-BB with history
prevfValuesForIstaBB = ql1_fValue(g,b,tau,x);
gPrevGlobal = g;
xPrevGlobal = x;
cgStatus='';

%% Outputs
xOut = x;
gOut = g;
bestOptimalityMeasure = optimalityMeasure(g,b,tau,x);

%% The actual algorithm loop
while 1
    %% Stopping conditions
    if optimalityMeasure(g,b,tau,x) <= accuracy
        xOut = x;
        gOut=g;
        algStatus='optimal';
        break;
    end
    
    if numA>=maxA
        algStatus='maxA';
        break;
    end
    
    %% First order step(s)
    if separate
        if gb(g,tau,x)
            %% First order step in subspace (x -> xF, g-> gF)
            gforfirstorderstep = g;
            gforfirstorderstep(x==0)=0;
            if numOuterIterations>0
                [xR, gR, numA,numB,fullHistory, prevfValuesForIstaBB,xPrevOutput] = ql1_istastep_bb(problem, gforfirstorderstep,x,alphabar,x-xPrevGlobal,g-gPrevGlobal,optimalityMeasure,accuracy,numA, numB,maxA, fullHistory, nargout, prevfValuesForIstaBB,M,2*xi,nu,xPrevOutput,outputLevel,stallingEpsilon);
            else
                xR=max(x-alphabar*(gforfirstorderstep+tau),0) - max(-x-alphabar*(-gforfirstorderstep+tau),0);
                gR=Ax(xR)-b;
                numA = numA+1;
                if (nargout>=2),fullHistory=alg_sub_RecordMV(fullHistory, numA,gF,b,tau,xR,optimalityMeasure,1);end
            end
            xPrevGlobal=x;
            gPrevGlobal=g;
        else
            if numOuterIterations>0
                [xR, gR, numA,numB,fullHistory, prevfValuesForIstaBB,xPrevOutput] = ql1_istastep_bb(problem, g,x,alphabar,x-xPrevGlobal,g-gPrevGlobal,optimalityMeasure,accuracy,numA,numB, maxA, fullHistory, nargout, prevfValuesForIstaBB,M,2*xi,nu,xPrevOutput,outputLevel,stallingEpsilon);
            else
                xR=max(x-alphabar*(g+tau),0) - max(-x-alphabar*(-g+tau),0);
                gR=Ax(xR)-b;
                numA = numA+1;
                if (nargout>=2),fullHistory=alg_sub_RecordMV(fullHistory, numA,gR,b,tau,xR,optimalityMeasure,1);end
            end
            xPrevGlobal=x;
            gPrevGlobal=g;
        end
    else
        %% Full space first order step  (x -> xR, g -> gR)
        if numOuterIterations>0
            [xR, gR, numA,numB,fullHistory, prevfValuesForIstaBB,xPrevOutput] = ql1_istastep_bb(problem, g,x,alphabar,x-xPrevGlobal,g-gPrevGlobal,optimalityMeasure,accuracy,numA,numB, maxA, fullHistory, nargout, prevfValuesForIstaBB,M,2*xi,nu,xPrevOutput,outputLevel,stallingEpsilon);
        else
            xR=max(x-alphabar*(g+tau),0) - max(-x-alphabar*(-g+tau),0);
            gR=Ax(xR)-b;
            numA = numA+1;
            if (nargout>=2),fullHistory=alg_sub_RecordMV(fullHistory, numA,gR,b,tau,xR,optimalityMeasure,1);end
        end
        xPrevGlobal=x;
        gPrevGlobal=g;
    end
    
    %% Output: after first order step computed
    if (outputLevel >=2),xPrevOutput = alg_sub_OutputX(gR,b,tau,xR,optimalityMeasure,numA,numB,'xR',(nargout >= 3), xPrevOutput) ; end
    
    %% Check for improvement
    if  optimalityMeasure(gR,b,tau,xR) < bestOptimalityMeasure
        xOut = xR;
        gOut = gR;
        bestOptimalityMeasure =optimalityMeasure(gOut,b,tau,xOut);
    end
    
    %% Stopping conditions
    if optimalityMeasure(gR,b,tau,xR) <= accuracy
        xOut = xR;
        gOut = gR;
        algStatus='optimal';
        break;
    end
    
    if numA>=maxA
        algStatus='maxA';
        break;
    end
    
    if (guessOptimal==1 && numOuterIterations>0 && isequal(sign(xR),workingOrthant) && (strcmp(cgStatus,'stall d') || strcmp(cgStatus,'stall dAd') || strcmp(cgStatus,'q inc')) )
        algStatus='guessOptimal';
        break;
    end
    
    %% CG (xE -> xCG, gE+tau.*workingOrthant->r)
    
    workingOrthant = sign(xR);
    xCG=xR;
    r=gR+tau.*workingOrthant;
    
    stepsSinceGoodCGpoint=0;
    cgNumMV=0;
    xG=xCG;
    rG=r ;
    prevguaranteedFminuscvvalue=ql1_fValue(r - tau.*workingOrthant,b,tau,xCG);
    currentQvalue=ql1_fValue(r - tau.*workingOrthant,b,tau,xCG);
    
    P = abs(workingOrthant);
    rho=P.*r;
    d=-rho;
    beta_CG=0;
    rrho=r'*rho;
    while 1
        if  ~gb(r- tau.*workingOrthant,tau,xCG)
            cgStatus = 'g balance';
            break;
        end
        d=-rho+beta_CG*d;
        if norm(d) <stallingEpsilon
            cgStatus='stall d';
            break;
        end
        Ad=Ax(d);
        cgNumMV=cgNumMV+1;
        numA=numA+1;
        if d'*Ad <stallingEpsilon
            if (nargout>=2),fullHistory=alg_sub_RecordMV(fullHistory, numA,Ax(xCG)-b,b,tau,xCG,optimalityMeasure,3);end
            cgStatus='stall dAd';
            break;
        end
        alpha=rrho/(d'*Ad);
        prevQvalue = currentQvalue;
        xCG=xCG+alpha*d;
        if (nargout>=2),fullHistory=alg_sub_RecordMV(fullHistory, numA,Ax(xCG)-b,b,tau,xCG,optimalityMeasure,3);end
        r=r+alpha*Ad;
        currentQvalue = 1/2 * xCG' * (r + b- tau.*workingOrthant) + (-b+tau.*workingOrthant)' * xCG;
        if (outputLevel >=3),xPrevOutput = alg_sub_OutputX(r - tau.*workingOrthant,b,tau,xCG,optimalityMeasure,numA,numB,'|xCG',(nargout >= 3),xPrevOutput) ; end
        if ((size(find(workingOrthant.*xCG<0),1)==0) ||(ql1_fValue(r - tau.*workingOrthant,b,tau,xCG) <= prevguaranteedFminuscvvalue  ))
            xG=xCG;
            rG=r;
            prevguaranteedFminuscvvalue=ql1_fValue(rG - tau.*workingOrthant,b,tau,xG) - c *norm(ql1_v(rG - tau.*workingOrthant,tau,xG))^2 ;
            xPrevGlobal=xCG-alpha*d;
            gPrevGlobal=r-alpha*Ad- tau.*workingOrthant;
            stepsSinceGoodCGpoint=0;
        else
            stepsSinceGoodCGpoint=stepsSinceGoodCGpoint+1;
        end
        if  optimalityMeasure(r - tau.*workingOrthant,b,tau,xCG) < bestOptimalityMeasure
            xOut = xCG;
            gOut= r - tau.*workingOrthant;
            bestOptimalityMeasure  =optimalityMeasure(r - tau.*workingOrthant,b,tau,xCG);
        end
        if optimalityMeasure(r - tau.*workingOrthant,b,tau,xCG) <= accuracy
            xOut = xCG;
            gOut = r - tau.*workingOrthant;
            cgStatus = 'optimal';
            break;
        end
        if  (  ql1_fValue(r - tau.*workingOrthant,b,tau,xCG) > prevguaranteedFminuscvvalue)
            cgStatus = 'F no dec';
            break;
        end
        if  (currentQvalue>prevQvalue)
            cgStatus = 'q inc';
            break;
        end
        rho=P.*r;
        beta_CG=r'*rho/(rrho);
        rrho=r'*rho;
    end
    
    %% Output: after CG procedure
    
    if (outputLevel >=2)
        fprintf('        --------cgNumMV = %i stepsSinceGoodCGpoint = %i cgStatus = %s\n', cgNumMV,stepsSinceGoodCGpoint, cgStatus);
        xPrevOutput = alg_sub_OutputX(rG - tau.*workingOrthant,b,tau,xG,optimalityMeasure,numA,numB,'xG',(nargout >= 3), xPrevOutput);
    end
    
    %% Record CG history
    
    numOuterIterations=numOuterIterations+1;
    
    if nargout>=2
        fullHistory.reasonForCGstop{numOuterIterations} = cgStatus;
        fullHistory.CGmvcount(numOuterIterations) = cgNumMV;
        fullHistory.CGglobalMVlink(numOuterIterations) = numA;
    end
    
    %% Check for improvement
    if  optimalityMeasure(rG - tau.*workingOrthant,b,tau,xG) < bestOptimalityMeasure
        xOut = xG;
        gOut=rG - tau.*workingOrthant;
        bestOptimalityMeasure =optimalityMeasure(gOut,b,tau,xOut);
    end
    
    %% Stopping conditions
    if  optimalityMeasure(rG - tau.*workingOrthant,b,tau,xG) <= accuracy
        xOut = xG;
        gOut=rG - tau.*workingOrthant;
        algStatus='optimal';
        break;
    end
    
    if rrho>stallingEpsilon && strcmp(cgStatus,'stall dAd') && (size(find(d.*xCG<0),1)==0)
        algStatus='unbounded';
        break;
    end
    
    if numA>=maxA
        algStatus='maxA';
        break;
    end
    
    %% CG post-processing (xCG -> xP,  r - tau.*workingOrthant->gP)
    
    % Checking that only a single step has been done in CG since the good
    % point, and that the good point was in the starting orthant. The stall
    % d check is used because in that case, Ad is not yet computed!
    if stepsSinceGoodCGpoint ==1 && ~strcmp(cgStatus,'stall d') && (size(find(xG.*workingOrthant<0),1)==0)
        ratioArray = xG./d;
        alphaf= max(ratioArray(ratioArray < 0));
        if size(alphaf,1)~=1 % Should not really happen, only extreme numerical problems cause this
            xP=xG;
            r=rG;
        else
            xP=xG-alphaf*d;
            xP(ratioArray==alphaf)=0;
            r=rG-alphaf*Ad;
            xPrevGlobal=xG;
            gPrevGlobal=rG- tau.*workingOrthant;
        end
    else
        xP=xG;
        r=rG;
    end
    
    gP = r - tau.*workingOrthant;
    
    %% Output: after post-processing
    if (outputLevel >=2),xPrevOutput = alg_sub_OutputX(gP,b,tau,xP,optimalityMeasure,numA,numB,'xP',(nargout >= 3),xPrevOutput) ; end
    
    %% Check for improvement
    if  optimalityMeasure(gP,b,tau,xP) < bestOptimalityMeasure
        xOut =xP;
        gOut=gP;
        bestOptimalityMeasure =optimalityMeasure(gOut,b,tau,xOut);
    end
    
    %% Check for stalling (if no progress is made througout the whole iteration)
    if norm(x-xP)<stallingEpsilon && norm(g-gP)<stallingEpsilon && size(find(x==0),1)==size(find(xP==0),1)
        algStatus='stall';
        break;
    end
    
    %%  (xP -> x,  gP->g)
    x=xP;
    g=gP;
    
end

%% Record and wrap up outputs
if outputLevel>=1
    xPrevOutput=alg_sub_OutputX(gOut,b,tau,xOut,optimalityMeasure,numA,numB,'xOut',(nargout >= 3), xPrevOutput);
    if numB==0
        fprintf('algStatus  = %s  numA  = %i  numCG  = %i\n',algStatus,numA,numOuterIterations);
    else
        fprintf('algStatus  = %s  numA  = %i  numB  = %i  numCG  = %i\n',algStatus,numA,numB,numOuterIterations);
    end
    fprintf('-----Finished alg_ql1-----  ');
    fprintf('\n');
end

varargout{1}.x = xOut;
varargout{1}.g = gOut;
varargout{1}.fValue = ql1_fValue(gOut,b,tau,xOut);
varargout{1}.optimalityMeasure = optimalityMeasure(gOut,b,tau,xOut);
varargout{1}.numZeros = size(find(xOut==0),1);
varargout{1}.numA = numA;
varargout{1}.numB = numB;
varargout{1}.numCGcycles = numOuterIterations;
varargout{1}.algStatus = algStatus;

if nargout >= 2
    fullHistory.fValues = fullHistory.fValues(1:numA);
    fullHistory.sparsity =  fullHistory.sparsity(1:numA);
    fullHistory.MVTypes =  fullHistory.MVTypes(1:numA);
    fullHistory.optimalityMeasures =  fullHistory.optimalityMeasures(1:numA);
    fullHistory.normV= fullHistory.normV(1:numA);
    fullHistory.kkterror= fullHistory.kkterror(1:numA);
    fullHistory.CGmvcount =  fullHistory.CGmvcount(1:numOuterIterations);
    fullHistory.CGglobalMVlink =  fullHistory.CGglobalMVlink(1:numOuterIterations);
    fullHistory.reasonForCGstop =  fullHistory.reasonForCGstop(1:numOuterIterations);
    varargout{2}=fullHistory;
end

if nargout >= 3
    varargout{3} = xPrevOutput;
end

end
