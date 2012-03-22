function solveHF(obj,envs,eps)
% envs: list of environments in which to solve HF (0=no environment)
% eps: tolerance for the HF convergence
if (nargin < 2)
   envs = 0:obj.nenv;
end
if (nargin < 3)
   eps = 1.0e-8;
end

nenv = size(envs,2);
for i = 1:nenv
   ienv = envs(i);
   if (ienv == 0)
      [obj.orb,obj.Eorb,obj.Ehf] = obj.hartreeFock(0,eps);
   else
    [obj.orbEnv(:,:,ienv), obj.EorbEnv(:,ienv), obj.EhfEnv(1,ienv)] = ...
       obj.hartreeFock(ienv,eps);
   end
end


end

