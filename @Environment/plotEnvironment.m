function plotEnvironment( obj, frag, rdisp )
%PLOTENVIRONMENT Visualize a fragment in a given environment.
%   Inputs:
%       obj: Environment object
%       frag: data from a single fragment, i.e. HL{1,1}
%       rdisp: displacement vector with data used to center the molecule
%            at (0,0), based on frag.rcart
%   Output:
%       3d Plot, where charges are shown as spheres proportional to the
%           magnitude of charge. Red corresponds to positive charges and
%           blue to negative charges. The molecule is shown for reference.

% rcart given in Bohr radii. Approx conversion factor.
rcartAng = frag.rcart / 1.889726124565062;

[x,y,z] = sphere(12);
%disp( sort( obj.rho ) );
%disp( sort( abs( obj.rho ) ) );
rhoMax = max( abs( obj.rho ) );

plot3( 0, 0, 0, '*' );
hold on;

% Draw the charges, proportional to their magnitude.
for icharge = 1:size( obj.r, 2 )
    scale = abs( obj.rho( icharge ) / rhoMax );
    a = x * scale;
    b = y * scale;
    c = z * scale;
    if obj.rho(icharge) >= 0
        color = zeros( size( x ) ) + 1;
    else
        color = zeros( size( x ) ) - 1;
    end
    surf( a + obj.r( 1, icharge ), b + obj.r( 2, icharge ), ...
        c + obj.r( 3, icharge ), color );
end

% For reference, draw the molecule as fixed-size circles.
for iatom = 1:size( rcartAng, 2 )
    plot3( rcartAng( 1, iatom ) + rdisp(1), ...
        rcartAng( 2, iatom ) + rdisp(2), ...
        rcartAng( 3, iatom ) + rdisp(3), 'o' );
end

%plot3( ctrChg(1), ctrChg(2), ctrChg(3), 'go' );

hold off;