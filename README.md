# Contents

This directory contains the Spectral Gradient Method with inexact projections described in:

O. P. Ferreira, M. Lemes, and L. F. Prudente. On the inexact scaled gradient projection method, technical report, 2021.

- File `SGMSDD.m` contains the function to solve 

   Minimize f(X)
   
   subject to X in { X in R^{nxn} : X = X^T, xii >= sum_{j~=i} |xij| for all i } 
   
   X in { X in R^{nxn} : xij >= 0 for all i,j }, 
   
   where f:R^{nxn} -> R is a continuous differentiable function. 
   
   The main algorithm uses the Dykstra's alternating projection algorithm to compute inexacts projections onto the feasible set, see file `Dykstra.m`.

- File `SGMSpec.m` contains the function to solve 
 
  Minimize f(X) 
  
  subject to X in { X in R^{nxn} : X = X^T, X >= 0, trace(X) = 1 }, 
  
  where f:R^{nxn} -> R is a continuous differentiable function. T
  
  The main algorithm uses the rank-k Frank-Wolfe algorithm to compute inexacts projections onto the feasible set, see file `SpecDp.m`.

- File `main.m` contains the main program, where you can choose the problem to be solved and the line search to be used. Modify files `evalf.m` and `evalg.m` to solve your own problem.

# Instructions:

The codes are written in Matlab. 

Execute the file `main.m` and see the output in the screen.

