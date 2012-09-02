function [ke, en, e2, newDensity] = modelCalcParallel(imod,ienv,updateDensity,...
   includeKE,includeEN,includeE2)
% Input:
%    imod          : int, model number
%    ienv          : int, environment number
%    updateDensity : bool, if density matrix needs to be recalculated
%    includeKE     : bool, calculate KE
%    includeEN     : vector of bool of 1:(maximum Z)

% retrieve model from mod
fileName = ['fitmeMod',num2str(imod),'.mat'];
load(fileName);
if (updateDensity)
   mod.solveHF(ienv);
   newDensity = mod.densitySave{ienv+1};
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
