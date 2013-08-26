function generateHLData(obj, spin)

    config = Fragment.defaultConfig();
    config.template = obj.HLtplFile;
    config.basisSet = obj.HLbasis;
    config.spin = spin;
    obj.HL = {};
    for igeom = 1:length(obj.geomPars)
       gpars = obj.geomPars{igeom};
       spars = ones(obj.HLnpar - length(gpars),1);
       config.par = [gpars(:); spars(:)];
       obj.HL{end+1} = Fragment(obj.dataDir, config, obj.cacheHL);
       for ienv = 1:length(obj.envs)
          obj.HL{end}.addEnv(obj.envs{ienv});
       end
    end
    obj.HLvalid = true;
end
