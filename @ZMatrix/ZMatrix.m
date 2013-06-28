classdef ZMatrix < handle
    %ZMatrix is a list of atoms from the ZAtom class
    %Pars can be a class that is constructed as needed.
    %Pars should be a list of bond lengths, angles and dihedral values
    %For more info on pars class, go to Pars_template class
    
    properties
        atoms = {};
        pars = {};      %Pars is in the following form, which does have to
                        %be it's own class. This could change if the 
                        %implementation is improved for auto-generating
                        %config files. Here is an example of what the
                        %class needs to have implemented for the 
                        %functions here to work:
                        %
                        %   pars.bond_pars = {}  List of bond vals in order
                        %   pars.ang_pars = {}   List of ang vals in order
                        %   pars.di_pars = {}    List of di vals in order
                        %
                        %Any other stuff the class has would be for the
                        %user, as the project needs.
                        
        file_name_base = '';
    end
    
    methods
        function zmat = ZMatrix( atoms, name )
            %Constructor
            %
            %file_name_base is automatically initialized to a random unique
            %name but can be changed as needed
            
            if nargin < 2
            	[~, name, ~] = fileparts( tempname() ); %Creates a random unique file name
                zmat.file_name_base = name;
            else
                zmat.file_name_base = name;
            end
            if nargin < 1
                zmat.atoms = {};
            else
                zmat.atoms = atoms;
            end
        end
        
        function make_atom( zmat, type, bond_ref, ang_ref, di_ref )
            %type should be a string
            %All inputs after type should be nonnegative integers
            %None of those inputs should be larger than the length of atoms
            %Has error messages for if conditions are not met (for range
            %and not having overlapping references, not for type checks)
            
            
            atom_num = length( zmat.atoms );
            if nargin > 3
                if di_ref == ang_ref
                    error('di_ref cannot equal ang_ref');
                end
                if di_ref == bond_ref
                    error('di_ref cannot equal ang_ref');
                end
                if di_ref > atom_num
                    error('di_ref must reference an existing atom in zmat.atoms');
                end
            else
                di_ref = 0;
            end
            
            if nargin > 2
                if ang_ref == bond_ref
                    error('ang_ref cannot equal bond_ref');
                end
                if ang_ref > atom_num
                    error('ang_ref must reference an existing atom in zmat.atoms');
                end
            else
                ang_ref = 0;
            end
            
            if nargin > 1
                if bond_ang > atom_num
                    error('bond_ref must reference an existing atom in zmat.atoms');
                end
            else
                bond_ref = 0;
            end
            
            if ~(nargin > 0)
                type = '';
            end
            
            atom = ZAtom( type, atom_num, bond_ref, ang_ref, di_ref );
            zmat.atoms{ atom_num } = atom;
        end
        
        function print_atoms( zmat, f_ID )
            %Calls print_atom for each atom in the zmat
            for i = 1:length( zmat.atoms )
                atom = zmat.atoms{i};
                atom.print_atom( f_ID );
            end
            fprintf( f_ID, '\n');
        end
        
        function print_pars_with_vars( zmat, f_ID)
            %Prints B, A, and D vals with the PAR instead of numbers
            num_atoms = length(zmat.atoms);
            num_bonds = num_atoms - 1;
            num_angs = num_atoms - 2;
            num_dis = num_atoms - 3;
            for i = 1:num_bonds
                fprintf( f_ID, ...
                    '   B%u             PAR%u\n', i, i);
            end
            for i = 1:num_angs
                fprintf( f_ID, ...
                    '   A%u             PAR%u\n', i, i + num_bonds);
            end
            for i = 1:num_dis
                fprintf( f_ID, ...
                    '   D%u             PAR%u\n', i, i + num_bonds + num_angs);
            end
            fprintf( f_ID, '\n!ENV' );
        end
                
        function print_z_to_tpl( zmat, path )
            full_file = fullfile( path, [ zmat.file_name_base, '.tpl' ] );
            f_ID = fopen( full_file, 'w' );
            %print_header( f_ID ); %THIS NEEDS TO BE DEFINED FOR YOUR PROJECT
            zmat.print_atoms( f_ID );
            zmat.print_pars_with_vars( f_ID );
            fclose( f_ID );
        end
        
        function print_config( zmat, path )
            %Creates a txt file to be used for Fragment
            
            name = strcat( zmat.file_name_base, '_config_raw' );
            full_file = fullfile( path, [ name, '.txt' ] );
            f_ID = fopen( full_file, 'w' );
            num_atoms = length(zmat.atoms);
            num_bonds = num_atoms - 1;
            num_angs = num_atoms - 2;
            num_dis = num_atoms - 3;
            for i = 1:num_bonds
                fprintf( f_ID, ...
                    '%d\n', zmat.pars.bond_pars{i} );
            end
            for i = 1:num_angs
                fprintf( f_ID, ...
                    '%d\n', zmat.pars.ang_pars{i});
            end
            for i = 1:num_dis
                fprintf( f_ID, ...
                    '%d\n', zmat.pars.di_pars{i});
            end
            fclose( f_ID );
        end
    end
    
end

