function solveHF(obj,eps)

if (nargin < 2)
   eps = 1.0e-8;
end

[obj.orb,obj.Eorb,obj.Ehf] = obj.hartreeFock(0,eps);

%
obj.EhfEnv  = zeros(1,obj.nenv);
obj.EorbEnv = zeros(obj.nbasis,obj.nenv);
obj.orbEnv  = zeros(obj.nbasis,obj.nbasis,obj.nenv);

for ienv = 1:obj.nenv
    [obj.orbEnv(:,:,ienv), obj.EorbEnv(:,ienv), obj.EhfEnv(1,ienv)] = ...
       obj.hartreeFock(ienv,eps);
end


end

