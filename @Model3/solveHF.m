function solveHF(obj,envs,eps,maxIter,minIter)
% envs: list of environments in which to solve HF (0=no environment)
% eps: tolerance for the HF convergence
if (nargin < 2)
   envs = 0:obj.nenv;
end
if (nargin < 3)
   eps = 1.0e-10;
end
if (nargin < 4)
   maxIter = 5000;
end
if (nargin < 5)
   minIter = 5;
end

nenv = size(envs,2);
for i = 1:nenv
   ienv = envs(i);
   if (ienv == 0)
      [obj.orb,obj.Eorb,obj.Ehf] = obj.hartreeFock(0,eps,maxIter,minIter);
   else
    [obj.orbEnv(:,:,ienv), obj.EorbEnv(:,ienv), obj.EhfEnv(1,ienv)] = ...
       obj.hartreeFock(ienv,eps,maxIter,minIter);
   end
end


end

