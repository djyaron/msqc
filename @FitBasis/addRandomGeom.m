function addRandomGeom(obj, ntimes)
if (nargin < 2)
   ntimes = 1;
end
npar = length(obj.geomRandomSpecs);
for ic=1:ntimes
   res = zeros(npar,1);
   for i = 1:npar
      t1 = obj.geomRandomSpecs{i};
      xlow = t1.low;
      xhigh = t1.high;
      res(i) = xlow + (xhigh-xlow).*rand(1,1);
   end
   obj.geomPars{end+1} = res;
end
end
