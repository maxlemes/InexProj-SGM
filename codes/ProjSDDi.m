function [x] = ProjSDDi(n,a,i)

% Given a vector a and index i, this function computes the projection of a
% onto the SSDi set, where
% SSDi = {A in S : aii >= sum_{j~=i}|aij|}.
%
% Reference:
% M. Raydan and P. Tarazaga, Primal and polar approach for computing the symmetric  
% diagonally dominant projection, Numer. Linear Algebra Appl. 9, pp.333-345, 2002.

% Define the element ai

ai = a(i);

% Define the set of indices ind that contains all indices from 1 to n except i

ind = setdiff([1:n],i);

% Define the vector aj, where aj = a with j ~= i

aj = a(ind);

% Test the conditions that guarantee that x=a or x=0 is the solution (see
% Theorem 3.1 of the above reference)

if ( ai >= sum( abs( aj ) ) )
    x = a;
    return
end

if ( ai < 0 && min( abs( ai ) > 2 * abs( aj ) ) == 1 )
    x = zeros(1,n);
    return
end

% Compute dbar (see Algorithm 3.2 of the above reference)

d = sum( abs( aj ) ) - ai;
c = n - nnz( ~aj ) + 1; % In the paper is missing the term +1
dbar = d/c;
s = 1;

while ( s == 1 )
    
    s = 0;
    
    for j = 1:n
        
        if ( j == i )
            continue
        end
        
        if ( a(j) > 0 )
            
            aux = a(j) - dbar;
            
            if ( aux < 0 )
                
                d = d - abs( a(j) );
                c = c - 1;
                a(j) = 0;
                s = 1;
                
            end
        end
        
        if ( a(j) < 0 )
            
            aux = a(j) + dbar;
            
            if ( aux > 0 )
                
                d = d - abs( a(j) );
                c = c - 1;
                a(j) = 0;
                s = 1;
                
            end
        end

    end
    
    dbar = d/c;
    
end

% Compute x
            
x(i) = ai + 2 * dbar;

for j = 1:n
    
    if ( j == i )
        continue
    end
    
    if ( a(j) == 0 )
        x(j) = 0;
    elseif ( a(j) > 0 )
        x(j) = a(j) - dbar;
    elseif ( a(j) < 0 )
        x(j) = a(j) + dbar;
    end
    
end