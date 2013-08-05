function res = Een(obj,iatom,envs)
% energy of interaction of the molecule with the environment
if (nargin < 3)
   envs = 0:obj.nenv;
end
n = size(envs,2);
res = zeros(1,n);
for i = 1:n
   ienv = envs(i);
   %res(i) = sum(sum( obj.density(ienv).*obj.H1en(:,:,iatom) ));
   res(i) = elementWiseCombine(obj.density(ienv), obj.H1en(:, :, iatom));
end
