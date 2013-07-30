function generateLLData(obj, spin)

    config = Fragment.defaultConfig();
    config.template = obj.LLtplFile;
    config.basisSet = 'GEN';
    config.spin = spin;
    obj.LL = {};
    spars = obj.extraPars;
    for igeom = 1:length(obj.geomPars)
       gpars = obj.geomPars{igeom};
       config.par = [gpars(:); spars(:)];
       obj.LL{end+1} = Fragment(obj.dataDir, config, obj.cacheLL);
       for ienv = 1:length(obj.envs)
          obj.LL{end}.addEnv(obj.envs{ienv});
       end
    end
    obj.LLvalid = true;
end