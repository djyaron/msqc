function fillInContexts(m1,envsTrain,m2,envsTest)
% m1 is the train set;  m2 is the test set
% load('datasets\ch4rDat.mat');
% m1 = LL(1:4,1);
% m2 = LL(5:8,1);
% envsTrain = 1:3;
% envsTest = 4:6;
%%
ntrain = length(m1);
% figure out how many atom types are in the train set
allTypes = [];
for i=1:ntrain
   allTypes = [allTypes,m1{i}.Z];
end
atypes = unique(allTypes);
%% determine the atom contexts
atomContexts = cell(length(atypes),1);
for itype = 1:length(atypes)
   atype = atypes(itype);
   ndim = sum(allTypes == atype) * length(envs);
   c1 = Context(ndim);
   % add all the data
   for imod = 1:length(m1)
      for ienv = envsTrain
         for iatom = find(m1{imod}.Z == atype)
            c1.addModel(m1{imod},ienv,iatom);
         end
      end
   end
   % do the fetaure extraction
   c1.extractFeatures;
   atomContexts{itype} = c1;
end

% determine the bond contexts

%% insert features into model based on the contexts
% want to save atomContextXSaved{iatom,ienv}
for imod = 1:length(m1)
   mm = m1{imod};
   mm.atomContextSaved = cell(mm.natom,mm.nenv + 1);
   for iatom = 1:m1{imod}.natom
      c1 = atomContexts{mm.aType(iatom)};
      for ienv = envsTrain
         mm.atomContextSaved{iatom,ienv+1} = c1.project(mm,iatom);
      end
   end
end
% do the same for the test set