function [X,f,time,outiter,NumEig,nfev,info] = SGMSpec(n,X,scale,LStype)

% This function is the Spectral Gradient Method with inexact projections
% to solve:
%
% Minimize f(X)
% subject to X in { X in R^{nxn} : X = X^T, X >= 0, trace(X) = 1 },
%
% where f:R^{nxn} -> R is a continuous differentiable function. The main 
% algorithm uses the rank-k Frank-Wolfe algorithm to compute inexacts 
% projections onto the feasible set.
%
% Reference:
% O. P. Ferreira, M. V. Lemes, and L. F. Prudente. On the inexact scaled 
% gradient projection method, technical report, 2021.
%
% Input variables:
%
% - n is the dimension of the problem
%
% - X is real nxn matrix
%    On entry X is the starting point
%
% - scale is a logical variable
%    scale = true:  the algorithm scales the objective function
%    scale = false: the algorithm does not use scalarization
%
% - LStype determines the line search to be used
%    LStype = 1: the algorithm uses the Armijo line search
%    LStype = 2: the algorithm uses the Average-type line search
%    LStype = 3: the algorithm uses the Max-type line search
%
% Output variables:
%
% - X is real nxn matrix
%    On exit X is the final iterate
%
% - f is the objective function value at the final iterate
%
% - time is the CPU time used
%
% - outiter is the number of outer iterations
%
% - Numeig is the number of computed eigenpairs
%
% - info contains the flag of the solution
%   info = 1: optimality satisfied
%   info = 2: the number of iterations was exhausted


global sF

% Start timing

tic;

% Initial Parameters

tol        = 10^(-6);
ftol       = 10^(-4);
maxoutiter = 1000;
alphamin   = 10^(-10);
alphamax   = 10^(10);
stpmin     = 10^(-10);
sigma1     = 0.1;
sigma2     = 0.9;
k          = 1;

gamma = 0.4999;

% Scale problem
	
scalefactor(n,X,scale);

% Counters

outiter = 0;
nfev    = 0;
FWit    = 0;
NumEig  = 0;

% Compute the function value

[f,flag] = sevalf(n,X);
nfev = nfev + 1;

if ( flag ~= 0 )
    info = -1;
    time = toc;
    return
end

% Compute the gradient

[G,flag] = sevalg(n,X);

if ( flag ~= 0 )
    info = -1;
    time = toc;
    return
end

normG = norm(G,'fro');

% Set the initial parameters for the line search

if ( LStype == 2 )
    q = 1;
    C = f;
    eta = 0.85;
end

if ( LStype == 3 )
    lastf(1) = f;
    M = 5;
end

% Define alpha

alpha = min(alphamax, max(alphamin, 1.0 / normG ) );

fprintf('--------------------------------------------------\n')
fprintf('  InexProj-SGM employing nonmonotone line search  \n')
fprintf('--------------------------------------------------\n')
fprintf('Dimension: %i x %i \n',n,n)
if ( LStype == 1 )
    fprintf('Line search: Armijo\n')
elseif ( LStype == 2 )
    fprintf('Line search: Average-type\n')
elseif ( LStype == 3 )
    fprintf('Line search: Max-type \n')
elseif ( LStype == 4 )
    fprintf('Line search: New nonmonotone line search \n')
end
fprintf('Tolerance for convergence: %.0e \n\n',tol)

while (1)  
    
    % Set the rank-k
    
    if ( outiter == 0 || outiter == 1 )
        easyProj = 0;
        k = 1;
    end
    
    if ( outiter > 2 )
        if ( ( infoProj == 0 || infoProj == 1 ) && outiterProj == 1 )
            easyProj = easyProj + 1;
            if ( easyProj == 3 )
                k = max( k - 1, 1 );
                easyProj = 0;
            end
        else
            easyProj = 0;
        end

    end
    
    % Project X - alpha * G onto the feasible set
    
    [Y,k,outiterProj,nCG,infoProj] = SpecDp(n,k,alpha,X,G,gamma);
    
    NumEig = NumEig + k + nCG;
    
    % Define the search direction
    
    d = Y - X;
    
    normd = max( max( abs( d ) ) );
    
    % Print information

    if ( mod(outiter,10) == 0 )
        fprintf('\n')
        fprintf('%-5s   %-8s   %2s     %-6s  %-3s    %-6s  %-8s %-8s\n','it','f','IS','nfev','k','nEig','||d||','|x-xprev|/|xprev|')
    end
    if ( outiter == 0 )
        fprintf('%5d   %5.2e   %2d   %6d   %3d   %6d   %8.2e        %1s\n',outiter,f,infoProj,nfev,k,NumEig,normd,'-')
    else
        fprintf('%5d   %5.2e   %2d   %6d   %3d   %6d   %8.2e     %8.2e\n',outiter,f,infoProj,nfev,k,NumEig,normd,normXXprev)
    end

    % --------------------------------
    % Stopping criteria
    % --------------------------------
    
    % Test whether the norm of d is too small
    
    if ( ( infoProj == 0 || infoProj == 1 ) && normd <= tol )
        info = 0;
        time = toc;
        
        % Print information
        
        fprintf('\n')
        fprintf('Solution was found.\n')
        fprintf('CPU time(s): %.1f \n',time)
        return
    end
    
    % Test whether the number of iterations is exhausted
    
    if ( outiter == maxoutiter )
        info = 1;
        time = toc;
        
        % Print information
        
        fprintf('\n')
        fprintf('The number of maximum iterations was reached.\n')
        fprintf('CPU time(s): %.1f \n',time)
        
        return
    end

    % --------------------------------
    
    % Increment outiter
    
    outiter = outiter + 1;
    
    % Define C
    
    if ( LStype == 1 )
        C = f;
    elseif ( LStype == 2 )
        if ( outiter == 1 )
            C = f;
        else
            q = eta * q + 1;
			C = ( q - 1 ) / q * C + f / q; 
        end
    elseif ( LStype == 3 )
        C = max( lastf );
    end
    
    % Compute the step size
    
    gtd = G(:)' * d(:);
    
    stp = 1.0;
    
    while (1)
        
        Xtrial = X + stp * d;
        
        [ftrial,flag] = sevalf(n,Xtrial);
        nfev = nfev + 1;
        
        if ( flag ~= 0 )
            info = -1;
            time = toc;
            return
        end

        if ( ftrial <= C + ftol * stp * gtd )
            break
        end
        
        if ( stp <= stpmin ) 
            break
            disp('Warning: stp = stpmin in the backtracking procedure')
        end
			
        stpq = ( ( gtd / ( ( f - ftrial) / stp + gtd ) ) / 2.0 ) * stp;

        if ( stpq >= sigma1 * stp && stpq <= sigma2 * stp )
            stp = stpq;
        else
            stp = stp / 2.0;
        end
        
    end
    
    % Update X
    
    Xprev = X;
    
    X = Xtrial;
    
    % Compute norm(X,Xprev)/norm(Xprev)
        
    normXXprev = norm( X - Xprev,'fro' )/norm( Xprev,'fro' );
    
    % Compute the function value

    f = ftrial;

    % Compute the gradient
    
    Gprev = G;

    [G,flag] = sevalg(n,X);
    
    if ( flag ~= 0 )
        info = -1;
        time = toc;
        return
    end

    normG = norm(G,'fro');
    
    % Store function value for the nonmonotone line search
    
    if ( LStype == 3 )
        lastf( mod( outiter, M ) + 1 ) = f;
    end

    % Define alpha

    s = X - Xprev;
    y = G - Gprev;
    
    a = s(:)' * s(:);
    b = s(:)' * y(:);
    
    if ( b <= 10^(-12) )
        alpha = alphamax;
    else
        alpha = min(alphamax, max(alphamin, a / b ) );
    end
        
end