function zmat = build_ethane_z( r1, r2, r3 )
    %creates a z matrix in the form C_H3_C_R1_R2_R3

    zmat = ZMatrix({});
    
    %Build static parts of ethane molecule
    zmat.make_atom( 6, 0, 0, 0 );
    zmat.make_atom( 6, 1, 0, 0 );
    zmat.make_atom( 1, 1, 2, 0 );
    zmat.make_atom( 1, 1, 2, 3 );
    zmat.make_atom( 1, 1, 2, 3 );
    
    %Build variable parts of ethane molecule
    add_group( zmat, r1 );
    add_group( zmat, r2 );
    add_group( zmat, r3 );
end

function add_group( zmat, r )
    %THIS FUNCTION ONLY WORKS FOR ETHANE IN LIMITED CASES
    %CAN ONLY BE USED FOR THE FORM C_H3_C_R1_R2_R3
    %ONLY WORKS WITH SINGLE ATOM R GROUPS AND LIMITED MULTIATOM R GROUPS

    bond = length( zmat.atoms );
    if is_cyano( r )
        zmat.make_atom( r{1}, 2, 1, 3 );
        zmat.make_atom( r{2}, bond + 1, 1, 3 );
    else
        for i = 1:length( r )
            if i == 1
                zmat.make_atom( r{i}, 2, 1, 3 );
                bond = bond + 1;
            else
                zmat.make_atom( r{i}, bond, 2, 1 );
            end
        end
    end
end

function bool = is_cyano( r )
   if length(r) ~= 2
       bool = false;
   else
       if  r{1} == 6 && r{2} == 7
           bool = true;
       else
           bool = false;
       end
   end
end