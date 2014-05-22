Ql1-Algorithm
=============

Optimization of L1 regularized quadratic 

A method for finding an optimal solution to
min (1/2)*x'*A*x - b'*x + norm(tau.*x,1)

Described here: http://www.ece.northwestern.edu/~nocedal/PDFfiles/ql1.pdf

Must have these inputs:

problem.Ax - function handle for matrix vector product

problem.b - vector b

problem.tau - vector or scalar tau



 Optional:
 problem.normA - An upper bound on the norm of matrix A
 (default = 1e6)
 
 opts.optimalityMeasure - Function handle for measuring optimality
 (default = scaled ISTA step)
 
 opts.accuracy - Accuracy that opt measure must reach
 (default = -1, see guessOptimal)
 
 opts.maxMV - Maximum number of Ax calls allowed
 (default = 10000)
 
 opts.outputLevel - How much output to display
      0 - none
      1 - Beginning and end
      2 - Outer loop progress (default)
      3 - Inner loop progress
      
 opts.x_0 - Starting point (default=all zeros)
 
 opts.separate - Separate into free vs zero steps
 0 no
 1 yes (default)
 
 opts.guessOptimal - Stop if guessing optimality
 0 no
 1 yes (default)





 Outputs:
 
 1st argument:
 
 x - Best solution found
 
 g - Gradient of the smooth part of F at x
 
 fValue - F evaluated at x
 
 optMeasure - Optimality measure at x
 
 numZeros - Sparsity of x
 
 numMV - Total work in term of Ax calls
 
 numOuterIterations - Number of outer iterations (CG cycles)
 
 algStatus - Final algorithm status. Can be:
        'optimal' - found solution to required accuracy
        'maxMV' - work limit reached
        'stall' - no progress during a full outer iteration
        'guessOptimal' - guess that opt solution found
 
 
 2nd argument: (optional and expensive! Only for debugging. Computed after each Ax call)
 
 fValues - F value of each iterate
 
 sparsity - Sparsity of each iterate
 
 MVTypes - What were Ax calls used for:
       0 - Initial g computation
       1 - Used for first order step
       2 - Used for relaxation step
       3 - Used for CG step
       
 optimalityMeasures - Optimality measures
 
 normV - Norms of minimum norm subgradient
 
 kkterror - KKT errors of reformulated problem
 
 CGmvcount - Ax calls for each CG cycle
 
 CGglobalMVlink - Current MV count after each CG cycle
 
 reasonForCGstop - Why CG cycles were stopped:
      'optimal' - Found optimal solution
      'g balance' - Gradient balance triggered
      'F no dec' - Not enough decrease in F
      'q inc' - Too much accuracy asked for
      'stall d' - Too much accuracy asked for
      'stall dAd' - Too much accuracy asked for
      
      
 3rd argument: (optional and expensive! Only for debugging. Computed after each Ax call)
      xPrevOutput - A sequence of all iterates
