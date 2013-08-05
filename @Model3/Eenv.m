function res = Eenv(obj,envs)
% energy of interaction of the molecule with the environment ienv
if (nargin < 2)
   envs = 0:obj.nenv;
end
n = size(envs,2);
res = zeros(1,n);
for i = 1:n
   ienv = envs(i);
   if (ienv == 0)
      res(i) = 0.0;
   else
      %res(i)  = sum(sum( obj.density(ienv).*obj.frag.H1Env(:,:,ienv) ) );
      res(i) = elementWiseCombine(obj.density(ienv), obj.frag.H1Env(:, :, ienv));
      res(i) = res(i) + obj.frag.HnucEnv(ienv) - obj.frag.Hnuc;
   end
end

