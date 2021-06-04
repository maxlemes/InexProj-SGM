function [X] = InitialPoint(n,Cset)

% Defining the initial point

if ( Cset == 1 )
    
    X = rand(n);
    X = ( X + X' )/2;

    for i = 1:n
        ind = setdiff([1:n],i);

        X(i,i) = 2 * sum( abs( X(i,ind) ) );
    end
    
    return
end

if ( Cset == 2 )
    
    X = 2 * rand(n) - 1;
    X = X * X';
    X = X / trace(X);
    
    return
end