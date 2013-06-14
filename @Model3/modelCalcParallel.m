function [ke, en, e2, newDensity, orb, Eorb, Ehf] = ...
   modelCalcParallel(obj, ienv, updateDensity, includeKE, includeEN, includeE2)
% Input:
%    imod          : int, model number
%    ienv          : int, environment number
%    updateDensity : bool, if density matrix needs to be recalculated
%    includeKE     : bool, calculate KE
%    includeEN     : vector of bool of 1:(maximum Z)

if (updateDensity)
   obj.solveHF(ienv);
   newDensity = obj.densitySave{ienv+1};
   if (ienv == 0)
      orb  = obj.orb;
      Eorb = obj.Eorb;
      Ehf  = obj.Ehf;
   else
      orb  = obj.orbEnv(:,:,ienv);
      Eorb = obj.EorbEnv(:,ienv);
      Ehf  = obj.EhfEnv(1,ienv);
   end
else
   newDensity = [];
end
if (includeKE == 1)
   ke = obj.EKE(ienv);
else
   ke = [];
end
en = zeros(obj.natom,1);
for iatom = 1:obj.natom
   if (includeEN(obj.Z(iatom)))
      en(iatom) = obj.Een(iatom,ienv);
   end
end
if (includeE2)
   e2 = obj.E2(ienv);
else
   e2 = 0.0;
end
