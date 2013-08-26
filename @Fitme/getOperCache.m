function res = getOperCache(obj)
nmod = length(obj.models);
res = cell(nmod,1);
for imod = 1:nmod
   res{imod} = obj.models{imod}.getCache;
end

end