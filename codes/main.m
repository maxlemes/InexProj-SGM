clear
clc

global A B c

scale = true;

seed = 123;
rng(seed);

% Define the constraint set
% Cset = 1: C = { X in R^{nxn} : X = X^T, xij >= 0 for all i,j, and xii >=
% sum_{j~=i} |xij| for all i }
% Cset = 2: C = { X in R^{nxn} : X = X^T, X >= 0, trace(X) = 1 }

Cset = 1;

% Define the line search to be used:
% LStype = 1: Armijo line search 
% LStype = 2: Average-type nonmonotone line search 
% LStype = 3: Max-type nonmonotone line search

LStype = 1;

% Define the dimension of the problem

n = 100;
m = 200; 

% Generating the objective function data

[A,B,c] = ObjFunction(n,m);

% Generating the starting point

[X0] = InitialPoint(n,Cset);

% Call the SGM solver

if ( Cset == 1 )
    [X,f,time,outiter,DIT,nfev,info] = SGMSDD(n,X0,scale,LStype);
end

if ( Cset == 2 )
    [X,f,time,outiter,DIT,nfev,info] = SGMSpec(n,X0,scale,LStype);
end