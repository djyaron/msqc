function [atypes, atomContexts, bondContexts] = ...
   fillInContexts(mtrain,envsTrain,mtest,envsTest,includeAdhoc)
% Input
%    mtrain      {ntrain} models for training
%    envsTrain   {ntrain}(1:nenv)  envs for each train model
%    mtest       {ntest}  models of test set
%    envsTest    {ntest}(1:nenv) envs for each test model
%    includeAdhoc   true to start context variables with Model3 adhoc
%                   (defaults to false)

if (nargin < 5)
   includeAdhoc = 0;
end
%%
ntrain = length(mtrain);
% determine atom types in the train set
allTypes = [];
for i=1:ntrain
   allTypes = [allTypes,mtrain{i}.aType];
end
atypes = unique(allTypes);

% determine the atom contexts
atomContexts = cell(length(atypes),1);
for itype = 1:length(atypes)
   atype = atypes(itype);
   % determine amount of training data
   ndim = 0;
   for imod=1:ntrain
      ndim = ndim + sum(mtrain{imod}.aType == atype) * ...
         length(envsTrain{imod});
   end
   c1 = Context(ndim,Context.atypeToZtype(atype));
   % add all the data
   for imod = 1:ntrain
      for ienv = envsTrain{imod}
         for iatom = find(mtrain{imod}.aType == atype)
            c1.addModel(mtrain{imod},ienv,iatom);
         end
      end
   end
   
   % do the fetaure extraction
   c1.extractFeatures;
   c1.plotLatent(1000+atype);
   atomContexts{itype} = c1;
end

%% verify project for atoms
ic = 0;
for imod = 1:ntrain
   for ienv = envsTrain{imod}
      for iatom = find(mtrain{imod}.aType == atype)
         ic = ic + 1;
         itype = find(atypes == mtrain{imod}.aType(iatom));
         c1 = atomContexts{itype};
         proj = c1.project(mtrain{imod},ienv,iatom);
         sc1 = c1.score(ic,:)';
         atomDiff(ic) = max(abs(proj - sc1));
      end
   end
end


%% determine the bond contexts
bondContexts = cell(length(atypes),length(atypes));
for type1 = 1:length(atypes)
   for type2 = type1:length(atypes)
      atype1 = atypes(type1);
      atype2 = atypes(type2);
      
      % determine amount of training data
      ndim = 0;
      for imod=1:ntrain
         atoms1 = find(mtrain{imod}.aType == atype1);
         atoms2 = find(mtrain{imod}.aType == atype2);
         for atom1 = atoms1
            for atom2 = atoms2
               if mtrain{imod}.isBonded(atom1,atom2)
                  ndim = ndim + length(envsTrain{imod});
               end
            end
         end
      end
      
      if (ndim == 0)
         continue
      end
      
      c1 = Context(ndim,Context.atypeToZtype(atype1), ...
         Context.atypeToZtype(atype2));
      % add all the data
      for imod = 1:ntrain
         for ienv = envsTrain{imod}
            for iatom = find(mtrain{imod}.aType == atype)
               atoms1 = find(mtrain{imod}.aType == atype1);
               atoms2 = find(mtrain{imod}.aType == atype2);
               for atom1 = atoms1
                  for atom2 = atoms2
                     if mtrain{imod}.isBonded(atom1,atom2)
                        c1.addModel(mtrain{imod},ienv,atom1,atom2);
                     end
                  end
               end
            end
         end
      end
      
      % do the fetaure extraction
      c1.extractFeatures;
      c1.plotLatent(2000 + atype1 + atype2);
      bondContexts{type1,type2} = c1;
      bondContexts{type2,type1} = c1;
      %end
   end
end

%% verify project for bonds
ic = 0;
for imod = 1:ntrain
   for ienv = envsTrain{imod}
      for iatom = find(mtrain{imod}.aType == atype)
         atoms1 = find(mtrain{imod}.aType == atype1);
         atoms2 = find(mtrain{imod}.aType == atype2);
         for atom1 = atoms1
            for atom2 = atoms2
               if (mtrain{imod}.isBonded(atom1,atom2))
                  ic = ic + 1;
                  itype = find(atypes == mtrain{imod}.aType(atom1));
                  jtype = find(atypes == mtrain{imod}.aType(atom2));
                  c1 = bondContexts{itype,jtype};
                  proj = c1.project(mtrain{imod},ienv,atom1,atom2);
                  sc1 = c1.score(ic,:)';
                  bondDiff(ic) = max(abs(proj - sc1));
               end
            end
         end
      end
   end
end

%% insert features into model based on the contexts
allMods = {mtrain{:},mtest{:}};
allEnvs = {envsTrain{:},envsTest{:}};
if (includeAdhoc)
   % get Model3 to fill in the ad hoc contexts, and then same them
   atomAdhoc = cell(length(allMods),1);
   bondAdhoc = cell(length(allMods),1);
   for imod = 1:length(allMods)
      mod = allMods{imod};
      mod.atomContext(1,1); % will fill in all atom contexts
      mod.bondContext(1,2,1); % will fill in all bond contexts
      atomAdhoc{imod} = mod.atomContextXSaved;
      bondAdhoc{imod} = mod.bondContextXSaved;
   end
end

for imod = 1:length(allMods)
   mod = allMods{imod};
   envs = allEnvs{imod};
   mod.atomContextXSaved = cell(mod.natom,mod.nenv + 1);
   for iatom = 1:mod.natom
      itype = find(atypes == mod.aType(iatom));
      c1 = atomContexts{itype};
      for ienv = allEnvs{imod}
         pcaContext = c1.project(mod,ienv,iatom);
         if (~includeAdhoc)
            mod.atomContextXSaved{iatom,ienv+1} = pcaContext;
         else
            adHoc = atomAdhoc{imod}{iatom,ienv+1};
            fullContext = [adHoc(:);pcaContext(:)];
            mod.atomContextXSaved{iatom,ienv+1} = fullContext;
         end
      end
   end
   
   mod.bondContextXSaved = cell(mod.natom,mod.natom,mod.nenv+1);
   for iatom = 1:mod.natom
      for jatom = find(mod.isBonded(iatom,:))
         itype = find(atypes == mod.aType(iatom));
         jtype = find(atypes == mod.aType(jatom));
         c1 = bondContexts{itype,jtype};
         for ienv = allEnvs{imod}
            pcaContext = c1.project(mod,ienv,iatom,jatom);
            if (~includeAdhoc)
               mod.bondContextXSaved{iatom,jatom,ienv+1} = ...
                  pcaContext;
            else
               adHoc = bondAdhoc{imod}{iatom,jatom,ienv+1};
               fullContext = [adHoc(:);pcaContext(:)];
               mod.bondContextXSaved{iatom,jatom,ienv+1} = ...
               fullContext;
            end
         end
      end
   end
end
