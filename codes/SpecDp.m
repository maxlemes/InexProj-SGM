function [W,k,outiter,NumEig,info] = SpecDp(n,k,alpha,X,G,gamma);

% Set up the vector to be projected

G = ( G + G' ) / 2;

Y = X - alpha * G;

% Counters

outiter = 0;
NumEig  = 0;

% Set initial w

W = X;

opts.issym  = 1;
opts.isreal = 1;
opts.disp   = 0;

%---------------------------------------------------------------------  
%     Main loop
%---------------------------------------------------------------------     

while (1)
    
    %---------------------------------------------------------------------  
    %     Stopping criteria
    %---------------------------------------------------------------------
    
    % Solve the CondG subproblem
    
    grad = Y - W;
    grad = ( grad + grad' ) / 2.0;
    grad = sparse( grad );
    
    opts.p = min( n, 60 );
    
    if ( isfield(opts,'v0') )
        opts = rmfield(opts,'v0');
    end

    [u,lambda,flagCG] = eigs( grad ,1,'la',opts);
    NumEig = NumEig + 1;
    
    opts = rmfield(opts,'p');
    
    if ( flagCG ~= 0 )
        disp('Warning: eigs failed to converge in the FW subproblem')
        gopt = -1e+20;
    end
    
    u( abs(u) < 1e-12 ) = 0;
    u = u / norm(u,2);

    Z = [];
    Z = u * u';

    Z = ( Z + Z' ) / 2.0;

    % Compute the optimal value of the subproblem

    gopt = real(sum(sum( (W - Y)' .* (Z - W) )));

    %  Test for convergence

    if ( outiter > 0 )
                
        phi = gamma * norm( W - X,'fro' )^2;
    
        if ( - gopt <= phi )
            info = 1;
            return
        end
    end
    
    if ( - gopt <= 10^(-15) )
        info = 0;
        return
    end
    
    if ( outiter > 0 && flagFW ~= 0 )
        info = 2;
        return
    end

    % Test whether the number of iterations is exhausted

    if ( k == n )
        info = 3;
        return
    end
    
    %-------------------------------------------------------------------
    %     Iteration
    %-------------------------------------------------------------------
    
    % Define k
    
    if ( outiter > 0 )  
        
        opts.v0 = randn(n,1);
        opts.v0(1:length(diagD)) = diagD;
        
        k = min( k + 1, n );
    end
    
    % Increment outiter

    outiter = outiter + 1;
    
    % Compute the k-Spectral Decomposition of Y
    
    opts.p = min( n, max( 3 * k , 50 ) );

    [U,D,flagFW] = eigs(Y,k,'la',opts);
    
    % Project the eigenvalues divided by eta in the simplex and save the
    % result in vector "sigma"
    
    diagD = diag(D);
    
    isNaNdiagD = find(isnan(diagD));
    
    diagD(isNaNdiagD) = [];
    U(:,isNaNdiagD) = [];
    
    if ( length( diagD ) == 0 )
        info = -2;
        return
    end
    
    k = length( diagD );

    [eigval,flag] = simplex_proj( diagD );   
    
    % Update W

    U( abs(U) < 1e-12 ) = 0.0;
    
    W = [];
    W = U * diag(eigval) * U';
    
    W = ( W + W' ) / 2.0;
    
%---------------------------------------------------------------------
%     Iterate
%---------------------------------------------------------------------    

end

%--------------------------------------------------------------------- 
%     End of main loop
%---------------------------------------------------------------------