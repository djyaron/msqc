function res = EKE(obj,envs)
% Kinetic energy in environment ienv
if (nargin < 2)
   envs = 1:obj.nenv;
end
n = size(envs,2);
res = zeros(1,n);
for ienv = 1:n
   res(ienv) = sum(sum( obj.density(ienv).*obj.KE(ienv) ) );
end

