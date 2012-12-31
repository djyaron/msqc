function [ke, en, e2, newDensity, orb, Eorb, Ehf] = modelCalcParallel(imod,ienv,updateDensity,...
   includeKE,includeEN,includeE2,scratchDir)
% Input:
%    imod          : int, model number
%    ienv          : int, environment number
%    updateDensity : bool, if density matrix needs to be recalculated
%    includeKE     : bool, calculate KE
%    includeEN     : vector of bool of 1:(maximum Z)

% retrieve model from mod
fileName = [scratchDir,'fitmeMod',num2str(imod),'.mat'];
load(fileName);
if (updateDensity)
   mod.solveHF(ienv);
   newDensity = mod.densitySave{ienv+1};
   if (ienv == 0)
      orb  = mod.orb;
      Eorb = mod.Eorb;
      Ehf  = mod.Ehf;
   else
      orb  = mod.orbEnv(:,:,ienv);
      Eorb = mod.EorbEnv(:,ienv);
      Ehf  = mod.EhfEnv(1,ienv);
   end
else
   newDensity = [];
end
if (includeKE == 1)
   ke = mod.EKE(ienv);
else
   ke = [];
end
en = zeros(mod.natom,1);
for iatom = 1:mod.natom
   if (includeEN( mod.Z(iatom) ))
      en(iatom) = mod.Een(iatom,ienv);
   end
end
if (includeE2)
   e2 = mod.E2(ienv);
else
   e2 = 0.0;
end
