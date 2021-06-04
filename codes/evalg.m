function [g] = evalg(n,X)

% Evaluating the objective function gradient

global A B c

g = A' * ( A * X - B );

x = diag(X);

for i = 1:n-1    
    g(i,i) = g(i,i) - 4 * c * x(i) * ( x(i+1) - x(i)^2 ) - 2 * ( 1 - x(i) );
    g(i+1,i+1) = g(i+1,i+1) + 2 * c * ( x(i+1) - x(i)^2 );
end