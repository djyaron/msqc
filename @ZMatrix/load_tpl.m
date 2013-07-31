function load_tpl( zmat, tpl_path, pars_path )
    %Takes in a tpl and pars txt and makes zmat have it's data
    %This is for backwards compatability

    if nargin < 1
        error( 'tpl path needed' );
    end
    
    tpl = fileread( tpl_path );
    atom_types = get_atom_types( tpl );
    bond_list = get_refs( tpl, 'B');
    ang_list = get_refs( tpl, 'A' );
    di_list = get_refs( tpl, 'D' );
    
    %Make atom list
%    zmat.atoms = cell(length(atom_types),1); Removed for temp fix.. needs
%    to have ZMatrix.make_atom updated to work with filling empty cells
    for i = 1:length(atom_types)
        if i == 1
            zmat.make_atom( atom_types{i} );
        elseif i == 2
            zmat.make_atom( atom_types{i}, bond_list(i-1) );
        elseif i == 3
            zmat.make_atom( atom_types{i}, bond_list(i-1), ...
                ang_list( i-2 ) );
        else
            zmat.make_atom( atom_types{i}, bond_list(i-1), ...
                ang_list( i-2 ), di_list(i-3) );
        end
    end
    
    %Make pars list if given
    if nargin < 3
        zmat.pars = {};
    else
        zmat.pars = get_pars( pars_path, length(atom_types));
    end
end

function atoms = get_atom_types( tpl )
    raw_atoms = regexp( tpl, '\s\w*ATOM', 'match' );
    atoms = cell( length(raw_atoms), 1 );
    for i = 1:length(raw_atoms)
        atom = char( raw_atoms(i) );
        atoms{i} = atom(2:length(atom)-4);
    end
end

function pars = get_pars( pars_path, n_atoms )
    f_ID = fopen( pars_path );
    pars_list = fscanf( f_ID, '%f' );
    fclose( f_ID );
    pars = {};
    bond_end = n_atoms - 1;
    ang_start = bond_end + 1;
    ang_end = bond_end + n_atoms - 2;
    di_start = ang_end + 1;
    di_end = ang_end + n_atoms - 3;
    pars.bond_pars = pars_list(1:bond_end);
    pars.ang_pars = pars_list(ang_start:ang_end);
    pars.di_pars = pars_list(di_start:di_end);
end

function list = get_refs( tpl, ref_type )
    
    exp = ['\s\d+\s', ref_type];
    raw_list = regexp( tpl, exp, 'match' );
    str_list = regexp( raw_list, '\d+', 'match' );
    list = zeros( length(str_list), 1 );
    for i = 1:length(str_list)
        list(i) = str2double( str_list{i} );
    end
end