function [e0,jacobian] = errJac(obj, p0)
% error and jacobian, calculated using cached operators
h = obj.epsDerivative;
obj.setPars(p0);
obj.cset.saveIndices;
obj.cacheAllOper;
e0 = obj.err(p0);
if (nargout > 1)
   save1 = obj.getOperCache;
   n = length(p0);
   jacobian = zeros(length(e0),n);
   for ip = 1:n
      p1 = p0;
      p1(ip) = p1(ip) + h;
      obj.setOperCache(save1);
      obj.setPars(p1);
      obj.cacheAllOper;
      e1 = obj.err(p1);
      jacobian(:,ip) = (e1-e0)/h;
   end
end
obj.clearOperCache;
   
end