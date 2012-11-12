function res = dipoleEnvironment( frag, cubeExtension, mag, ndipole )
%DIPOLEENVIRONMENT Double cube approach to environment creation.
%   Inputs:
%       frag: data from a single fragment, i.e. HL{1,1}
%       cubeExtension [ x y z ]: amount to extend box beyond the inner box
%           that bounds the molecule
%       mag: randomly generated charges are in range -magnitude..magnitude
%       ndipole: the number of charges to place in the environment
%   Output:
%       res: Environment object

% Data to be used in determining minimum safe distance for charges.
load( [ '@Environment', filesep, 'atomicRadii.mat' ] );

rcartAng = frag.rcart; 

res = Environment;
icharge = 1;
ncharge = ndipole * 2;

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
badDipoles = 0;
while icharge <= ncharge
    % Randomly place a point charge.
    testLoc1 = [ 0; 0; 0 ];
    testLoc1( 1 ) = rand() * outerBox( 1, randi( [ 1 2 ] ) );
    testLoc1( 2 ) = rand() * outerBox( 1, randi( [ 1 2 ] ) );
    testLoc1( 3 ) = rand() * outerBox( 1, randi( [ 1 2 ] ) );
    % Create a dipole relative to first charge in polar coordinates.
    testDip = [ 0; 0; 0 ]; % [ r; theta; phi ]
    testDip( 1 ) = rand();
    testDip( 2 ) = rand() * pi;
    testDip( 3 ) = rand() * 2 * pi;
    % Convert to absolute Cartesian coordinates.
    testLoc2 = [ 0; 0; 0 ];
    testLoc2( 1 ) = testLoc1( 1 ) + testDip( 1 ) * sin( testDip( 2 ) ) ...
        * cos( testDip( 3 ) );
    testLoc2( 2 ) = testLoc1( 2 ) + testDip( 1 ) * sin( testDip( 2 ) ) ...
        * sin( testDip( 3 ) );
    testLoc2( 3 ) = testLoc1( 3 ) + testDip( 1 ) * cos( testDip( 2 ) );
    
    chgGood = true;
    for iatom = 1:frag.natom
        dist1 = sqrt( sum( ( testLoc1 - rcartAng( :, iatom ) ) .^ 2 ) );
        dist2 = sqrt( sum( ( testLoc2 - rcartAng( :, iatom ) ) .^ 2 ) );
        % Minimun safe distances based on staying 5A from H and 6A from C.
        safeDist = 7.143 * atomicRadii( frag.Z( iatom ) ) + 1.214;
        if safeDist > dist1 || safeDist > dist2
            chgGood = false;
            break;
        end
    end
    if chgGood == true
        r( :, icharge ) = testLoc1;
        r( :, icharge + 1 ) = testLoc2;
        icharge = icharge + 2;
    else
        badDipoles = badDipoles + 1;
    end
end

% Generate charge magnitudes.
rhoParity = randi( [ 1 2 ], 1, ndipole );
rhoParity = rhoParity - mod( rhoParity, 2 ) - 1;
rhoScale = rand( 1, ndipole );
rho = rhoParity .* rhoScale .* mag;
ncharge = 2 * ndipole;
res.rho = zeros( 1, ncharge );
res.rho( 1:2:ncharge ) = rho;
res.rho( 2:2:ncharge ) = rho;

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
    disp( badDipoles );

    disp( 'Center of Charge: [ x; y; z ]' );

    disp( ctrChg );

    disp( 'Inside inner box? [ x; y; z ]' );

    disp( bound );

    disp( 'Magnitude of Electric Field at Origin' );

    disp( transpose( sum( transpose( OField ) ) ) );

    res.plotEnvironment( frag, [ 0 0 -0.72 ], ctrChg );

end



end
