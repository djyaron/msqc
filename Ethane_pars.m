classdef Ethane_pars < handle
    %Hacky class that gives values to the pars of ethane
    %
    
    properties
        %Where the par values are stored in the form of key value pairs
        bond_pars = {};
        ang_pars = {};
        di_pars = {};
        
        %Estimate values that are used to fill pars
        bond_vals = {'C' 'C' 1.54
                     'C' 'H' 1.14
                     'C' 'O' 1.43
                     'C' 'N' 1.47
                     'C' 'F' 1.35
                     'C' 'Cl' 1.77
                     'N' 'H' 1.001
                     'N' 'O' 1.25
                     'O' 'H' 0.958
                    };
        ang_vals = {120.00001   %Bent
                    120.00001   %Tri
                    109.5 %Tetrahedral
                   };
        di_vals = {0.0001, 120.00001, -120.00001, 180.000001};
    end
    
    methods
        function eth = Ethane_pars( atoms )
            %Constructor
            eth.get_bond_vals( atoms );
            eth.get_ang_vals( atoms );
            eth.get_di_vals( atoms );
        end
        
        function val = get_bond_val( eth, a_type1, a_type2 )
            %Finds the value of the bond between the 2 atoms
            %Linear searches the atom vals for the matching line
            for i = 1:length( eth.bond_vals )
                if a_type1 == eth.bond_vals{i, 1}
                    if a_type2 == eth.bond_vals{i, 2}
                        val = eth.bond_vals{i, 3};
                        break
                    end
                elseif a_type1 == eth.bond_vals{i, 2}
                    if a_type2 == eth.bond_vals{i, 1}
                        val = eth.bond_vals{i, 3};
                        break
                    end
                end
            end
        end
        
        function get_bond_vals( eth, atoms )
            %Starts at 2 because that is the first atom with a bond_ref
            for i = 2:length(atoms)
                atom1 = atoms{i};
                atom2 = atoms{atom1.bond_ref};
                eth.bond_pars{i-1} = eth.get_bond_val( atom1.type, atom2.type );
            end
        end
        
        function val = get_ang_val( eth, atoms, atom_i )
            bonded_ref = atoms{ atom_i }.bond_ref;
            num_bonded = atoms{ bonded_ref }.bond_total;
%             for i = atom_bonded:length(atoms)
%                 if atoms{i}.bond_ref == atom_bonded
%                    num_bonded = num_bonded + 1; 
%                 end
%             end
            val = eth.ang_vals{num_bonded - 1};
        end
        
        function get_ang_vals( eth, atoms )
            for i = 3:length(atoms)
                eth.ang_pars{i-2} = eth.get_ang_val( atoms, i );
            end
        end
        
        function val = get_di_val( eth, atoms, atom_i )
            %di_ref = atoms{ atom_i }.di_ref;
            if atom_i == 4
                val = eth.di_vals{2};
            elseif atom_i == 5
                val = eth.di_vals{3};
            elseif atom_i == 6
                val = eth.di_vals{1};
            elseif atoms{atom_i}.bond_ref == 2
                bool = false;
                for i = atom_i+1:length(atoms)
                    if atoms{i}.bond_ref == 2
                        bool = true;
                        break
                    end
                end
                if bool
                    val = eth.di_vals{2};
                else
                    val = eth.di_vals{3};
                end
            elseif atoms{atom_i-1}.bond_ref == 2
                val = 0;
            else
                val = 180;
            end
        end
        
        function get_di_vals( eth, atoms )
            for i = 4:length(atoms)
                eth.di_pars{i-3} = eth.get_di_val( atoms, i );
            end
        end
        
        function print_file( eth, file_base )
            count = 1;
            name = strcat( file_base, '_config.txt' );
            f_ID = fopen( name, 'w' );
            
            for i = 1:length(eth.bond_pars)
                fprintf( f_ID, 'PAR%u\t%f\n', count, eth.bond_pars{i} );
                count = count + 1;
            end
            for i = 1:length(eth.ang_pars)
                fprintf( f_ID, 'PAR%u\t%f\n', count, eth.ang_pars{i} );
                count = count + 1;
            end
            for i = 1:length(eth.di_pars)
                fprintf( f_ID, 'PAR%u\t%d\n', count, eth.di_pars{i} );
                count = count + 1;
            end
            fclose( f_ID );
        end
    end
    
end

