function [x,flag] = simplex_proj(a)

% function [x,flag] = simplex_proj(a)
%
% Projects a vector 'a' onto the unit simplex

n = length(a);

flag = 1;

x = a;

I = 1:n;

for t = 1:n
    d = length(I);
    
    x(I) = x(I) + (1/d) * (1 - sum(x(I)));
    
    N = find( x < 0 );

    if ( isempty(N) )
        flag = 0;
        return
    end

    x(N) = 0;
    
    I = find( x > 0 );
end