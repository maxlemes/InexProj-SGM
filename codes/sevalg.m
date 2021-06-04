function [G,flag] = sevalg(n,X) 

global sF

flag = 0;

[G] = evalg(n,X);

for i = 1:n
    for j = 1:n
        [isnum] = IsANumber( G(i,j) );

        if ( isnum == 0 )
            flag = -1;
            disp('WARNING: There is an element whose value may be +Inf, -Inf or NaN')
            disp('in the gradient of the objective computed by the user-supplied')
            disp('subroutine evalg.')
        end
    end
end

G = sF * G;