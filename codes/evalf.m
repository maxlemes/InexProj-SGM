function [f] = evalf(n,X)

% Evaluating the objective function:
% f(X) = 1/2 * |AX-B|_F^2 + sum_{i=1}^{n-1}[c(X(i+1,i+1)-X(i,i)^2)^2+(1-X(i,i))^2]

global A B c

f = 0.5 * norm( A * X - B,'fro')^2;

x = diag(X);

for i = 1:n-1
    f = f + c * ( x(i+1) - x(i)^2 )^2 + ( 1 - x(i) )^2;
end