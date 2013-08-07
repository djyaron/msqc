clear( 'all' );

path = 'D:\small_test_data';
r_groups = [
    1
    9
    ];

%Loops through each r group for the 3 r groups
for i = 1:length(r_groups)
    r1 = r_groups(i);
    for j = 1:i
        r2 = r_groups(j);
        for k = 1:j
            r3 = r_groups(k);
            for l = 1:k
                disp( [i, j, k, l] );
                r4 = r_groups(l);
                
                zmat = ZMatrix();
                zmat.make_atom( 6, 0, 0, 0 );
                zmat.make_atom( r1, 1, 0, 0 );
                zmat.make_atom( r2, 1, 2, 0 );
                zmat.make_atom( r3, 1, 2, 3 );
                zmat.make_atom( r4, 1, 2, 3 );
                
                zmat.file_name_base = 'methane';
                zmat.pars.bond_pars = {1.2, 1.2, 1.2, 1.2};
                zmat.pars.ang_pars = {109.5, 109.5, 109.5};
                zmat.pars.di_pars = {120.00001, -120.00001};
                
                temp_path = tempname([path,'\fluorine_testdat']);
                mkdir( temp_path );
                
                config = Config();
                config.opt = 1;
                config.zmat = zmat;
                
                frag = Fragment();
                bool = frag.runFragment( temp_path, config );
                
                if bool == 1
                    vary_mol( temp_path, frag );
                end
            end
        end
    end
end
