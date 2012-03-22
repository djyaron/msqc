function res = EKE(obj,envs)
% Kinetic energy in environment ienv
if (nargin < 2)
   envs = 0:obj.nenv;
end
n = size(envs,2);
res = zeros(1,n);
for i = 1:n
   ienv = envs(i);
   res(i) = sum(sum( obj.density(ienv).*obj.KE(ienv) ) );
end

