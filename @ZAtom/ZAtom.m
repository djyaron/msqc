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
        
        function print_atom( atom, f_ID )
            %Outputs line of z matrix for atom
            %Prints according to if there is a bond, ang or dihedral
            
            str = strcat(atom.type);
            fprintf( f_ID, ' %sATOM%u', str, atom.num);
            if atom.num > 1
                fprintf( f_ID, ' %u B%u', atom.bond_ref, atom.num - 1);
            end
            if atom.num > 2
                fprintf( f_ID, '    %u A%u', atom.ang_ref, atom.num - 2);
            end
            if atom.num > 3
                fprintf( f_ID, '    %u D%u    0', atom.di_ref, atom.num - 3);
            end
            fprintf( f_ID, '\n');
        end
        
        function up_bond_total( atom )
            atom.bond_total = atom.bond_total + 1;
        end
    end
    
end

