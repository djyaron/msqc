function vary_mol( path )

    zip = ls( [path, '\*.zip'] );
    if isempty( zip )
        return
    end
    
    tpl = [path, '\ethane.tpl'];
    base_dir = [path, '\vary'];
    max_iter = 5;
    

    unzip( [path, '\', zip], path );
    old_config = [path, '\full_opt_config.txt'];
    
    for i = 1:max_iter
        new_dir = [base_dir, num2str( i )];
        if exist( new_dir, 'file' ) ~= 7
            mkdir( new_dir );
            new_config = [new_dir, '\config_vary_', num2str( i ), '.txt' ];
            copyfile( tpl, new_dir );
            pars = rand_config( old_config, new_config );
            %run_gaus( new_dir, pars, 'STO-3G' );
            run_gaus( new_dir, pars, '631-G' );
        end
    end
end


function run_gaus( path, pars, basisSet )

    config = Fragment.defaultConfig();
    config.template = 'ethane';
    config.par = pars;
    config.basisSet = basisSet;
    Fragment( path, config );
end

function varied_nums = rand_config( old_config, new_config )
    %Varies config file for a molecule from the optimized config

    config_ID = fopen( old_config );
    nums = fscanf( config_ID, '%f' );
    fclose( config_ID );
    start_ang_i = start_angs( nums );
    varied_nums = zeros( length(nums), 1 );
    
    for i = 1:length(nums)
        temp = floor( i/start_ang_i );
        varied_nums(i) = vary_val( nums(i), temp);
    end
    new_ID = fopen( new_config, 'w' );
    for i = 1:length(varied_nums)
        fprintf( new_ID, '%f\n', varied_nums(i) );
    end
    fclose( new_ID );
end

function new = vary_val( num, type )
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

function pos = start_angs( num_list )
    %Gives position of angles in list assuming the first angle has to be
    %larger than bond length values
    pos = 0;
    max_bond = 10;
    for i = 1:length( num_list )
        if num_list(i) > max_bond
            pos = i;
            break
        end
    end
end