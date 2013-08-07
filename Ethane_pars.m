classdef Ethane_pars < handle
    %Hacky class that gives values to the pars of ethane
    %
    
    properties
        %Where the par values are stored in the form of key value pairs
        bond_pars = {};
        ang_pars = {};
        di_pars = {};
        
        %Estimate values that are used to fill pars
        %bond_vals is in the form:
                    %Z Z bond_length
        bond_vals = {6 6 1.54
                     6 1 1.14
                     6 8 1.43
                     6 7 1.47
                     6 9 1.35
                     6 17 1.77
                     7 1 1.001
                     7 8 1.25
                     8 1 0.958
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
                eth.bond_pars{i-1} = eth.get_bond_val( atom1.z, atom2.z );
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
        
        function vary_pars( eth )
            for i = 1:length(eth.bond_pars)
                eth.bond_pars{i} = eth.vary_val( eth.bond_pars{i}, 0);
            end
            for i = 1:length(eth.ang_pars)
                eth.ang_pars{i} = eth.vary_val( eth.ang_pars{i}, 1);
            end
            for i = 1:length(eth.di_pars)
                eth.di_pars{i} = eth.vary_val( eth.di_pars{i}, 1);
            end
        end
        
        function new = vary_val( eth, num, type )
            %Outputs varied value for bond or angle
            %num is value to be varied
            %type must be 0 for bond and 1 for angle
            
            new = (-1)^randi(2) * rand(); %random number between [-1, 1]
            bond_vary = 0.3;
            ang_vary = 10;
            
            if type
                new = num + new * ang_vary;
            else
                new = num + new * bond_vary;
            end
        end
    end
    
end

