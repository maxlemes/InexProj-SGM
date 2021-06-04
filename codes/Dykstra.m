function [X,alphamax,a,it,info] = Dykstra(n,lambda,X0,G,eta)

% This function is the Dykstra's alternating projection algorithm to obtain 
% the projection onto the primal cone SDD+. The algorithm projects cyclically
% onto the n convex sets SDD+i whose non-empty intersection is the cone SDD+.
%
% References:
% [1] E. G. Birgin, J. M. Martinez, and M. Raydan, Inexact spectral projected 
% gradient methods on convex sets. IMA J. Numer. Anal., 23(4):539?559, 2003.
%
% [2] M. Raydan and P. Tarazaga, Primal and polar approach for computing the symmetric  
% diagonally dominant projection, Numer. Linear Algebra Appl. 9, pp.333-345, 2002.
%
% See Algorithm 3.1 of [1] and Algorithm 4.1 of [2]

% Parameters:

tol  = 10^(-8);

stopit = 2;

% Set up the matrix to be project onto SDD+

G = ( G + G' ) / 2;

A = X0 - lambda * G;

% Initialize the increments

sp = issparse(A);

if ( sp )
    z = sparse(n,n);
    zbox = sparse(n,n);
else
    z = zeros(n);
    zbox = zeros(n);
end

% Counters:

it = 0;

% Initialize variable c, a, Q, and the auxiliar vector dotzy

c  = 0;
a0 = - lambda^2 * norm(G,'fro')^2;
a  = a0;

Q = - Inf;

dotzy(1:n+1) = 0;

% Compute the auxiliar vector aux

for i = 1:n
    ind = setdiff([1:n],i);
    
    xiimxind(i) = X0(i,i) - sum( abs( X0(i,ind) ) );
end

%----------------------------------
% Main algorithm    
%----------------------------------
    
while (1)
    
    % Test convergence
    
    if ( - a <= tol )

        info = 0;
        
        X = ( X + X' ) / 2;
        return
    end
    
    if ( it > 0 )
        
        % Compute Q(X-X0)
        
        Qprev = Q;
        
        Q = norm( X - X0 , 'fro' )^2 / ( 2 * lambda ) + sum( sum( G .* ( X - X0 ) ) );
    
        if ( Q <= eta * a )
            info = 1;            
            X = ( X + X' ) / 2;
            return
        end
    end
    
    if ( it > 0 && abs( Q - Qprev )/abs(Q) <= 10^(-6) && norm( A - Aprev2 ,'fro' ) <= tol  )
        stop = stop + 1;
    else
        stop = 0;
    end
    
    if ( it > 0 && stop == stopit )
        info = 2;
        X = ( X + X' ) / 2;
        return
    end
    
    % Iterate
    
    it = it + 1;
    
    Aprev2 = A;
    
    %----------------------------------
    % Beginning of a Dykstra iteration
    %----------------------------------
    
    for i = 1:n
        
        % Save the i-th line of A in y

        y = A(i,:);
        
        % Compute yz as y minus the increment
                
        yz = y - z(i,:);
        
        % Projects yz onto SSDi+
        
        [u] = ProjSDDi(n, full(yz),i);
        
        if ( sp )
            u = sparse(u);
        end
        
        % Update matrix A
        
        A(i,:) = u;
        A(:,i) = u;
        
        % Incorporate the current information into c
        
        ind = setdiff([1:n],i);
        
        c = c + ( y(i) - u(i) )^2 + 2 * sum( ( y(ind) - u(ind) ).^2 );
        c = c + 2 * ( z(i,i) * u(i) + 2 * sum( z(i,ind) .* u(ind) ) - dotzy(i) );
        
        % Update the increment
        
        z(i,:) = A(i,:) - yz;
        
        % Compute <z_i,y_i> to be used in the next iteration
        
        dotzy(i) = z(i,i) * u(i) + 2 * sum( z(i,ind) .* u(ind) );
    end
    
    % Projecting onto the box constraint
    
    Aprev = A;
    
    A = max( A - zbox, 0 );
    
    % Incorporate the current information into c
    
    c = c + norm( A - Aprev,'fro' )^2;
    c = c + 2 * ( sum( sum( zbox .* A ) ) - dotzy(n+1) );
    
    % Update the increment
    
    zbox = A - ( Aprev - zbox );
    
    % Compute <z_i,y_i> to be used in the next iteration
    
    dotzy(n+1) = sum( sum( zbox .* A ) );
  
    %----------------------------------
    % End of a Dykstra iteration
    %----------------------------------
    
    % Update a
    
    a = ( c + a0 ) / ( 2 * lambda );
    
    % Compute alphamax
    
    D = A - X0;
    
    alphamax = 10^(20);
    
    for i = 1:n
        
        ind = find( D(i,:) < - 10^(-10) );
        if ( length(ind) > 0 )
            alphamax = min( alphamax, min( - X0(i,ind) ./ D(i,ind) ) );
        end
        
        ind = setdiff([1:n],i);
        
        diimdind = D(i,i) - sum( D(i,ind) );
        
        if ( diimdind < - 10^(-10) )
            alphamax = min( alphamax, - xiimxind(i) / diimdind );
        end
    end
    
    % Update X
    
    X = X0 + min( alphamax, 1 ) * D;

    %----------------------------------
    % End of main algorithm    
    %----------------------------------
    
end