function scalefactor(n,X,scale)

global sF

eps = 1.0d-08;

% Compute objective function scaling factor

if ( scale )

    [G] = evalg(n,X);
    
    sF = 1 / max( 1, norm( G, 'inf' ) );
    sF = max( eps, sF );

else
    sF = 1;
end