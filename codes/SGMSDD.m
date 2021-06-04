function [X,f,time,outiter,DIT,nfev,info] = SGMSDD(n,X,scale,LStype)

% This function is the Spectral Gradient Method with inexact projections
% to solve:
%
% Minimize f(X)
% subject to X in { X in R^{nxn} : X = X^T, xii >= sum_{j~=i} |xij| for all i }
%            X in { X in R^{nxn} : xij >= 0 for all i,j },
%
% where f:R^{nxn} -> R is a continuous differentiable function. The main 
% algorithm uses the Dykstra's alternating projection algorithm to compute 
% inexacts projections onto the feasible set.
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
% - DIT is the number of Dykstra?s iterations
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

zeta = 0.80;
beta = 0.85;

% Scale problem
	
scalefactor(n,X,scale);

% Counters

outiter = 0;
nfev    = 0;
DIT     = 0;

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
    
    % Project X - alpha * G onto the feasible set

    [Y,alphaStpmax,a,init,infoProj] = Dykstra(n,alpha,X,G,zeta);
        
    DIT = DIT + init;
    
    normd = max( max( abs( Y - X ) ) );

    if ( infoProj ~= 0 )

        % Define the search direction

        if ( alphaStpmax >= 1 / beta )
            d = Y - X;
        elseif ( alphaStpmax >= 1 )
            d = beta * alphaStpmax * ( Y - X );
        else
            d = beta * ( Y - X );
        end 
        
    end
    
    % Print information

    if ( mod(outiter,10) == 0 )
        fprintf('\n')
        fprintf('%-5s   %-8s   %2s    %3s   %-6s %-8s     %-2s    %-8s\n','it','f','IS','DIT','nfev','||d||','-a','|x-xprev|/|xprev|')
    end
    if ( outiter == 0 )
        fprintf('%5d   %5.2e   %2d %6d %6d  %8.2e   %8.2e       %1s\n',outiter,f,infoProj,DIT,nfev,normd,-a,'-')
    else
        fprintf('%5d   %5.2e   %2d %6d %6d  %8.2e   %8.2e    %8.2e\n',outiter,f,infoProj,DIT,nfev,normd,-a,normXXprev)
    end

    % --------------------------------
    % Stopping criteria
    % --------------------------------
    
    % Test whether the norm of d is too small
    
    if ( infoProj ~= 2 && normd <= tol )
        info = 1;
        time = toc;
        
        % Print information
        
        fprintf('\n')
        fprintf('Solution was found.\n')
        fprintf('CPU time(s): %.1f \n',time)
        return
    end
    
    % Test whether the Dykstra algorithm has indentified the optimality of X
    
    if ( infoProj == 0 )
        info = 1;
        time = toc;
        
        % Print information
        
        fprintf('\n')
        fprintf('Solution was found.\n')
        fprintf('CPU time(s): %.1f \n',time)
        return
    end
      
    % Test whether the number of iterations is exhausted
    
    if ( outiter == maxoutiter )
        info = 2;
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