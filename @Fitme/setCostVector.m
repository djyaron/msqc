function setCostVector(obj,eps)

if (nargin < 2)
   eps = 0.1;
end
% cache derivative with respect to every parameter in the mixers,
% independent of whether they are current fixed. 
if (isempty(obj.mixerCost))
   obj.mixerCost = initializeMixerCost(obj,eps);
end

% Pull out the cost for nonfixed value within the mixers
res = zeros(1,obj.npar);
ic = 1;
for imix = 1:length(obj.mixers)
   mix = obj.mixers{imix};
   np = mix.npar;
   if (np > 0)
      ifree = find(mix.fixed == 0);
      res(ic:(ic+np-1)) = obj.mixerCost{imix}(ifree);
      ic = ic + np;
   end
end

obj.costVector = res/norm(res);

end

function mixerCost = initializeMixerCost(obj,eps)
% Don't want to include cost at this point
saveCost = obj.cost;
obj.cost = 0.0;
nmix = length(obj.mixers);
% save all params, and set temporarily to zero
psave = cell(nmix,1);
for imix = 1:nmix
   mix = obj.mixers{imix};
   psave{imix} = mix.par;
   mix.par = zeros(size(mix.par));
end
% get derivative with respect to every parameter
err0 = obj.err(obj.getPars);
mixerCost = cell(nmix,1);
for imix = 1:nmix
   mix = obj.mixers{imix};
   npar = length(mix.par);
   t1 = zeros(1,npar);
   for ip = 2:npar
      mix.par(ip) = 0.1;
      err1 = obj.err(obj.getPars);
      t1(ip) = norm(err1-err0)/eps;
      mix.par(ip) = 0.0;
   end
   mixerCost{imix} = abs(t1);
end
% restore 
obj.cost = saveCost;
for imix = 1:nmix
   mix = obj.mixers{imix};
   mix.par = psave{imix};
end

end
