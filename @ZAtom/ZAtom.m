classdef ZAtom < handle
    %REPRESENTS AN ATOM IN A Z MATRIX
    %PROPERTIES INCLUDES ALL INTERNAL REFERENCE INFO
	%INFO STORED AS POINTERS TO OTHER ATOMS IN THE Z MATRIX LIST
    %PROPERTY FOR NUMBER OF OTHER ATOMS BONDED TO THIS ATOM
    %METHOD TO OUTPUT TO A LINE OF A Z MATRIX FILE
    
    properties
        type = ''; %carbon ('C'), hydrogen ('H'), etc
        num = 0; %absolute reference number for atom in molecule
        
        bond_ref = 0; %Zatom.num of atom to create bond from 
                      %self to bond_ref
                      %0 if does not exist
        
        ang_ref = 0; %Zatom.num of atom to create angle from
                     %self to bond_ref to ang_ref
                     %0 if does not exist
        
        di_ref = 0; %Zatom.num of atom to create dihedral from
                    %self, bond_ref, ang_ref and di_ref
                    %0 if does not exist;
                    
        bond_total = 0; %number of atoms bonded to this atom
    end
    
    methods
        function atom = ZAtom( type, num, bond_ref, ang_ref, di_ref)
            %Constructor method]
            %FOR USAGE USE ARGS AS NEEDED, DO NOT FILL 0s FOR DEFAULTS
            %Uses 0 for ref values if not applicable
            
            if nargin < 5
                atom.di_ref = 0;
            else
                atom.di_ref = di_ref;
            end
            if nargin < 4
                atom.ang_ref = 0;
            else
                atom.ang_ref = ang_ref;
            end
            if nargin < 3
                atom.bond_ref = 0;
                atom.bond_total = 0;
            else
                atom.bond_ref = bond_ref;
                if bond_ref ~= 0
                    atom.bond_total = 1;
                else
                    atom.bond_total = 0;
                end
            end
            if nargin < 2
                atom.num = 0;
            else
                atom.num = num;
            end
            if nargin < 1
                atom.type = '';
            else
                atom.type = type; 
            end
        end
        
        function text = atom_text( atom, bq )
            %Outputs line of z matrix for atom
            %Prints according to if there is a bond, ang or dihedral
            
            text = '';
            str = strcat(atom.type);
            text = [text, ' ', str];
            
            if nargin > 1
                if bq == 1
                    text = [text, '-Bq'];
                end
            end
            if atom.num > 1
                bond_str = num2str( atom.bond_ref );
                var_num_str = num2str( atom.num - 1 );
                text = [text, ' ', bond_str, ' B', var_num_str];
            end
            if atom.num > 2
                ang_str = num2str( atom.ang_ref );
                var_num_str = num2str( atom.num - 2 );
                text = [text, '    ', ang_str, ' A', var_num_str];
            end
            if atom.num > 3
                di_str = num2str( atom.di_ref );
                var_num_str = num2str( atom.num - 3 );
                text = [text, '    ', di_str, ' D', var_num_str];
            end
            text = [text, '\n'];
        end
        
        function up_bond_total( atom )
            atom.bond_total = atom.bond_total + 1;
        end
        
        function find_charge( atom )
            
        end
    end
    
end

