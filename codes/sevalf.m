function [f,flag] = sevalf(n,X) 

global sF

flag = 0;

[f] = evalf(n,X);

[isnum] = IsANumber(f);

if ( isnum == 0 )
    flag = -1;
    disp('WARNING: The objective function value computed by the')
    disp('user-supplied subroutine evalf may be +Inf, -Inf or NaN.')
end

f = sF * f;