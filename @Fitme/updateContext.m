function updateContext(obj)

for imod = 1:length(obj.models)
   mod = obj.models{imod};
   for ienv = obj.envs{imod}
      mod.updateContext(ienv);
   end
end

for imod = 1:length(obj.models)
   obj.cset.addModel( obj.models{imod} );
end

end
