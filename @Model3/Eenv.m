function res = Eenv(obj,envs)
% energy of interaction of the molecule with the environment ienv
if (nargin < 2)
   envs = 1:obj.nenv;
end
n = size(envs,2);
res = zeros(1,n);
for ienv = 1:n
   res(ienv)  = sum(sum( obj.density(ienv).*obj.frag.H1Env(:,:,ienv) ) );
   res(ienv) = res(ienv) + obj.frag.HnucEnv(ienv);
end

