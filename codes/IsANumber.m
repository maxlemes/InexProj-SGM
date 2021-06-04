function [isnum] = IsANumber(x)
  
bignum = 1.0d+99;

isnum = true;
if ( abs( x ) > bignum || isnan( x ) || isinf( x ) ) 
    isnum = false;
end