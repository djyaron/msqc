clear( 'all' );

path = 'C:\Users\Zack\Documents\GitHub\msqc\ethaneRs';
r_groups = {
    {'H'}
    {'F'}
    {'Cl'}
    {'O' 'H'}
    {'N' 'H' 'H'}
    {'N' 'O' 'O'}
    {'C' 'N'}
    {'C' 'O' 'H'}
    };

%Loops through each r group for the 3 r groups
for i = 1:length(r_groups)
    r1 = r_groups{i};
    for j = 1:length(r_groups)
        r2 = r_groups{j};
        for k = 1:length(r_groups)
            r3 = r_groups{k};
            zmat = build_ethane_z(r1, r2, r3);
            zmat.file_name_base = 'ethane';
            zmat.pars = Ethane_pars( zmat.atoms );
            temp_path = tempname([path,'\','testdat']);
            mkdir( temp_path );
            zmat.print_z_to_tpl( temp_path );
            zmat.print_config( temp_path );
        end
    end
end
