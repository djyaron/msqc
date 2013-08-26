function clearOperCache(obj)
for imod = 1:length(obj.models)
   for ienv = obj.envs{imod}
      obj.models{imod}.clearCache(ienv);
   end
end

end