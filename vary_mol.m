function vary_mol( path, frag )

    max_iter = 5;
    absolute_max = 8;
    base_zmat = frag.config.zmat;
    
    i = 1;
    while i < max_iter
        config = Config();
        config.basisSet = '6-31G';
        config.zmat = base_zmat;
        vary_pars( config.zmat.pars );
        
        frag = Fragment();
        bool = frag.runFragment( path, config );
        
        if bool == 1 && max_iter <= absolute_max
            max_iter = max_iter + 1;
        end
        i = i + 1;
    end
end

function vary_pars( pars )
    for i = 1:length(pars.bond_pars)
        pars.bond_pars{i} = vary_val( pars.bond_pars{i}, 0);
    end
    for i = 1:length(pars.ang_pars)
        pars.ang_pars{i} = vary_val( pars.ang_pars{i}, 1);
    end
    for i = 1:length(pars.di_pars)
        pars.di_pars{i} = vary_val( pars.di_pars{i}, 1);
    end
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