function res = Een(obj,iatom,envs)
% energy of interaction of the molecule with the environment
if (nargin < 3)
   envs = 1:obj.nenv;
end
n = size(envs,2);
res = zeros(1,n);
for ienv = 1:n
   res(ienv) = sum(sum( obj.density(ienv).*obj.H1en(iatom,ienv) ));
end
