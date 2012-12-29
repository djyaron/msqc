function updateContext(obj)

for imod = 1:length(obj.models)
   mod = obj.models{imod};
   % do no env, since we are taking charge differences
   mod.updateContext(0);
   for ienv = obj.envs{imod}
      mod.updateContext(ienv);
   end
end

for imod = 1:length(obj.models)
   obj.cset.addModel( obj.models{imod} );
end

end
