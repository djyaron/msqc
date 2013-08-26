function cacheAllOper(obj)

for imod = 1:length(obj.models)
   for ienv = obj.envs{imod}
      obj.models{imod}.cacheOperators(ienv);
   end
end

end