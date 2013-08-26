function res = plotBasis(obj,ifunc,fnum,rrange)
% ifunc = which basis function
% fnum = figure number (0 to not plot)
% rrange = range of r
   res = zeros(size(rrange));
   prim = obj.basisPrims{ifunc};
   for ip = 1:obj.basisNprims(ifunc)
      coef = prim(1,ip);
      expp = prim(2,ip);
      res = res + coef * exp(-expp * rrange.^2);
   end
   if (fnum > 0)
      figure(fnum)
      plot(rrange,res);
   end
end

