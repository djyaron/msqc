function [fitme,cset] = makeFitme(obj,ms,f1)

if (nargin < 3)
   HLfromF1 = false;
else
   HLfromF1 = true;
end

notIn = setdiff(ms.atomTypes,obj.atomTypes);
if (~isempty(notIn))
   error(['MFactory.makeFitme: factory created without these types',...
      num2str(notIn)]);
end

% Clear current, and then add all, modifiers
model = ms.models;
for imod = 1:length(model)
   model{imod}.clearModifiers;
end

for imix = 1:length(obj.mixInfo)
   minfo = obj.mixInfo{imix};
   for imod = 1:length(model)
      mod = model{imod};
      switch minfo.type
         case 'KEdiags'
            mod.addKEmodDiag(minfo.iatom,1,minfo.mixer);
         case 'KEdiagp'
            mod.addKEmodDiag(minfo.iatom,2,minfo.mixer);
         case 'ENdiags'
            mod.addENmodDiag(minfo.iatom,1,minfo.mixer);
         case 'ENdiagp'
            mod.addENmodDiag(minfo.iatom,2,minfo.mixer);
         case 'E2diag'
               mod.addH2modDiag(minfo.iatom,minfo.mixer);
         case 'KEbondss'
            mod.addKEmodBonded(minfo.iatom,minfo.jatom,1,1,minfo.mixer);
         case 'KEbondsp'
            mod.addKEmodBonded(minfo.iatom,minfo.jatom,1,2,minfo.mixer);
         case 'KEbondps'
            mod.addKEmodBonded(minfo.iatom,minfo.jatom,2,1,minfo.mixer);
         case 'KEbondpp'
            mod.addKEmodBonded(minfo.iatom,minfo.jatom,2,2,minfo.mixer);
         case 'KEbondh'
            mod.addKEmodBondedh(minfo.iatom,minfo.jatom,minfo.mixer);
         case 'ENbondss'
            mod.addENmodBonded(minfo.iatom,minfo.jatom,1,1,minfo.mixer);
         case 'ENbondsp'
            mod.addENmodBonded(minfo.iatom,minfo.jatom,1,2,minfo.mixer);
         case 'ENbondps'
            mod.addENmodBonded(minfo.iatom,minfo.jatom,2,1,minfo.mixer);
         case 'ENbondpp'
            mod.addENmodBonded(minfo.iatom,minfo.jatom,2,2,minfo.mixer);
         case 'ENbondh'
            mod.addENmodBonded1h(minfo.iatom,minfo.jatom,minfo.mixer);
         case 'E2bond'
            mod.addH2modOffDiag(minfo.iatom,minfo.jatom,minfo.mixer);
         otherwise
            error('Mfactory: unrecognized mixInfo type');
      end
   end
end

% Create a cset to hold all needed contexts
cset = CSet;
for imod = 1:length(model)
   cset.addModel(model{imod});
end
cset.saveIndices;
% Create the fitme object
fitme = Fitme;
for imod = 1:length(ms.models)
   if (HLfromF1)
      fitme.addFrag(ms.models{imod},[],ms.pnum(imod));
   else
      fitme.addFrag(ms.models{imod},ms.HLfrag{imod},ms.pnum(imod));
   end
end
fitme.includeKE = 1;
fitme.includeEN = ones(1,20);
fitme.includeE2 = 1;
fitme.includeEtot = 1;
fitme.silent = 0;
fitme.plot = 0;
fitme.parallel = 1;
if (HLfromF1)
   fitme.envs = f1.envs;
   fitme.HLKE = f1.HLKE;
   fitme.HLEN = f1.HLEN;
   fitme.HLE2 = f1.HLE2;
   fitme.LLKE = f1.LLKE;
   fitme.LLEN = f1.LLEN;
   fitme.LLE2 = f1.LLE2;
   fitme.parHF = [];
   fitme.HLs = [];
else
   fitme.setEnvs(ms.envs);
   % setEnvs calculates the HL values of everything we are fitting to, so the
   % HLs are no longer needed. By removing these from fitme, we make the fitme
   % object quite a bit smaller.
   fitme.HLs = [];
end
end

