function vary_mol( path, frag )

    max_iter = 5;
    base_zmat = frag.config.zmat;
    
    for i = 1:max_iter
        config = Config();
        config.basisSet = '6-31G';
        config.zmat = base_zmat;
        config.zmat.pars.vary_pars();
        
        frag = Fragment();
        frag.runFragment( path, config );
    end
end
