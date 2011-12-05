function res = newBox( frag, cubeExtension, mag, ncharge )
%NEWBOX Double cube approach to environment creation.
%   Inputs:
%       frag: data from a single fragment, i.e. HL{1,1}
%       cubeExtension [ x y z ]: amount to extend box beyond the inner box
%           that bounds the molecule
%       mag: randomly generated charges are in range -magnitude..magnitude
%       ncharge: the number of charges to place in the environment
%   Output:
%       res: Environment object

% Data to be used in determining minimum safe distance for charges.
load( [ '@Environment', filesep, 'atomicRadii.mat' ] );

% rcart given in Bohr radii. Approx conversion factor.
% Should really be changed in readfchk
rcartAng = frag.rcart / 1.889726124565062;

res = Environment;
icharge = 1;

% Find bounding box of the molecule.
innerBox = [ min( rcartAng( 1, : ) ) max( rcartAng( 1, : ) ); ...
             min( rcartAng( 2, : ) ) max( rcartAng( 2, : ) ); ...
             min( rcartAng( 3, : ) ) max( rcartAng( 3, : ) ); ];       
         
% Compute limits of outer box, based on inner box.
outerBox = zeros( 3, 2 );
outerBox( :, 1 ) = -1 * cubeExtension;
outerBox( :, 2 ) = cubeExtension;
outerBox = outerBox + innerBox;
         
rng('shuffle');         

r = zeros( 3, ncharge );

% Place charges until enough acceptable charges have been found.
% Currently, there is no check on 'acceptable.'
badCharges = 0;
while icharge <= ncharge
    testLoc = [ 0; 0; 0 ];
    testLoc( 1 ) = rand() * outerBox( 1, randi( [ 1 2 ] ) );
    testLoc( 2 ) = rand() * outerBox( 1, randi( [ 1 2 ] ) );
    testLoc( 3 ) = rand() * outerBox( 1, randi( [ 1 2 ] ) );
    chgGood = true;
    for iatom = 1:frag.natom
        dist = sqrt( sum( ( testLoc - rcartAng( :, iatom ) ) .^ 2 ) );
        % Minimun safe distances based on staying 5A from H and 6A from C.
        safeDist = 7.143 * atomicRadii( frag.Z( iatom ) ) + 1.214;
        if safeDist > dist
            chgGood = false;
            break;
        end
    end
    if chgGood == true
        r( :, icharge ) = testLoc;
        icharge = icharge + 1;
    else
        badCharges = badCharges + 1;
    end
end

% Generate charge magnitudes.
rhoParity = randi( [ 1 2 ], 1, ncharge );
rhoParity = rhoParity - mod( rhoParity, 2 ) - 1;
rhoScale = rand( 1, ncharge );
res.rho = rhoParity .* rhoScale .* mag;

res.ncharge = ncharge;
res.r = r;


if 0
    % Debugging:
    ctrChg = [ sum( r( 1, : ) .* res.rho ) / sum( abs( res.rho ) ); ...
               sum( r( 2, : ) .* res.rho ) / sum( abs( res.rho ) ); ...
               sum( r( 3, : ) .* res.rho ) / sum( abs( res.rho ) ) ]; 
       
    bound = ( ctrChg > innerBox( :, 1 ) ) == ( ctrChg < innerBox( :, 2 ) );

    dist = sqrt( r(1,:).^2 + r(2,:).^2 + r(3,:).^2 );
    e0 = 8.854197e-12;
    fieldMag = ( 1 / ( 4 * pi * e0 ) ) * res.rho .* ( 1.6e-19 ) ...
        ./ ( dist .^ 3 );
    OField = repmat( fieldMag, 3, 1 ) .* res.r;
    
    disp( 'Inner Box: [ min max ]' );        
    disp( innerBox );  

    disp( 'Outer Box: [ min max ]' )
    disp( outerBox );

    disp( '# Bad Charges:' );
    disp( badCharges );

    disp( 'Center of Charge: [ x; y; z ]' );

    disp( ctrChg );

    disp( 'Inside inner box? [ x; y; z ]' );

    disp( bound );

    disp( 'Magnitude of Electric Field at Origin' );

    disp( transpose( sum( transpose( OField ) ) ) );

    res.plotEnvironment( frag, [ 0 0 -0.72 ], ctrChg );

end



end

