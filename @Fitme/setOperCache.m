function setOperCache(obj,c1)
nmod = length(obj.models);
for imod = 1:nmod
   obj.models{imod}.setCache(c1{imod});
end

end