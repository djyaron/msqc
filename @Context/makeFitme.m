function [ftrain ftest] = makeFitme(mtrain,envsTrain,HLTrain,mtest, ...
   envsTest,HLTest,includeAdhoc,separateSP,include1s)
% Input
%    mtrain      {ntrain} models for training
%    envsTrain   {ntrain}(1:nenv)  envs for each train model
%    mtest       {ntest}  models of test set
%    envsTest    {ntest}(1:nenv) envs for each test model

if (nargin < 7)
   includeAdhoc = 0;
end
if (nargin < 8)
   separateSP = 0;
end
if (nargin < 9)
   include1s = 0;
end

% get all the necessary contexts for this set of molecules
[atypes, atomContexts, bondContexts] =  ...
   Context.fillInContexts(mtrain,envsTrain,mtest,envsTest,includeAdhoc);
%% create a mixer for every atom and bond context

% need to expand # context variables by 3 if including adhoc
if (includeAdhoc)
   extraContexts = 3;
else
   extraContexts = 0;
end

ic = 0;
mixInfo = cell(0,0);
for itype =1:length(atypes)
   atype = atypes(itype);
   if (~isempty(atomContexts{itype}))
      %Mixer(parIn,mixType,desc,funcType)
      ncontexts = atomContexts{itype}.ndim + extraContexts;
      
      desc = ['KE atype ',num2str(atype),' pca '];
      if (Context.atypeToZtype(atype) == 6)
         if (separateSP)
            functype = 2; % scale without constant
            parIn = [1 zeros(1,ncontexts)];
            fixed = [0 ones(1,ncontexts)]; % fix all contexts
            mixType = 11; % atom context mixer
            m1 = Mixer(parIn, mixType, desc, functype);
            m1.fixed = fixed; % fix all contexts
            minfo.mixer = m1;
            minfo.type = 'KEdiags';
            minfo.atype1 = atype;
            mixInfo{end+1} = minfo;
            minfo.type = 'KEdiagp';
            minfo.mixer = m1.deepCopy;
            mixInfo{end+1} = minfo;
         else
            functype = 3; % scale with constant
            parIn = [1 zeros(1,ncontexts) 0];
            fixed = [0 ones(1,ncontexts) 0]; % fix all contexts
            mixType = 11; % atom context mixer
            m1 = Mixer(parIn, mixType, desc, functype);
            m1.fixed = fixed; % fix all contexts
            minfo.mixer = m1;
            minfo.type = 'KEdiag';
            minfo.atype1 = atype;
            mixInfo{end+1} = minfo;
         end
      else
         functype = 2; % scale without constant
         parIn = [1 zeros(1,ncontexts)];
         fixed = [0 ones(1,ncontexts)]; % fix all contexts
         mixType = 11; % atom context mixer
         m1 = Mixer(parIn, mixType, desc, functype);
         m1.fixed = fixed; % fix all contexts
         minfo.mixer = m1;
         minfo.type = 'KEdiag';
         minfo.atype1 = atype;
         mixInfo{end+1} = minfo;
      end

      if (include1s && (Context.atypeToZtype(atype) == 6))
         desc = ['KE 1s1s ',num2str(atype),' pca '];
         mixType = 11; % atom context mixer
         functype = 2; % scale without constant
         parIn = [1 zeros(1,ncontexts)];
         m1 = Mixer(parIn, mixType, desc, functype);
         m1.fixed = [0 ones(1,ncontexts)]; % fix all contexts
         minfo.mixer = m1;
         minfo.type = 'KE1s1s';
         minfo.atype1 = atype;
         mixInfo{end+1} = minfo;
         
         desc = ['KE 1s2s ',num2str(atype),' pca '];
         mixType = 11; % atom context mixer
         functype = 2; % scale without constant
         parIn = [1 zeros(1,ncontexts)];
         m1 = Mixer(parIn, mixType, desc, functype);
         m1.fixed = [0 ones(1,ncontexts)]; % fix all contexts
         minfo.mixer = m1;
         minfo.type = 'KE1s2s';
         minfo.atype1 = atype;
         mixInfo{end+1} = minfo;
      end
      
      desc = ['EN atype ',num2str(atype),' pca '];
      if (Context.atypeToZtype(atype) == 6)
         functype = 3; % scale with constant
         parIn = [1 zeros(1,ncontexts) 0];
         fixed = [0 ones(1,ncontexts) 0]; % fix all contexts
      else
         functype = 2; % scale without constant
         parIn = [1 zeros(1,ncontexts)];
         fixed = [0 ones(1,ncontexts)]; % fix all contexts
      end
      mixType = 11; % atom context mixer
      m1 = Mixer(parIn, mixType, desc, functype);
      m1.fixed = fixed;
      minfo.mixer = m1;
      minfo.type = 'ENdiag';
      minfo.atype1 = atype;
      mixInfo{end+1} = minfo;
      
      desc = ['E2 atype ',num2str(atype),' pca '];
      mixType = 11; %
      functype = 2; % scale without constant
      parIn = [1 zeros(1,ncontexts)];
      fixed = [0 ones(1,ncontexts)]; % fix all contexts
      m1 = Mixer(parIn, mixType, desc, functype);
      m1.fixed = fixed;
      minfo.mixer = m1;
      minfo.type = 'E2diag';
      minfo.atype1 = atype;
      mixInfo{end+1} = minfo;
      for jtype = itype:length(atypes)
         atype2 = atypes(jtype);
         if (~isempty(bondContexts{itype,jtype}))
            ncontexts = bondContexts{itype,jtype}.ndim + extraContexts;
            parIn = [1 zeros(1,ncontexts)];
            mixType = 12; % off-diagonal context mixer
            desc = ['KE atypes ',num2str(atype),' ',num2str(atype2),' pca '];
            functype = 2; % scale without const
            fixed = [0 ones(1,ncontexts)]; % fix all contexts
            m1 = Mixer(parIn, mixType, desc, functype);
            m1.fixed = fixed;
            m1.hybrid = 1;
            minfo.mixer = m1;
            minfo.type = 'KEbond';
            minfo.atype1 = atype;
            minfo.atype2 = atype2;
            mixInfo{end+1} = minfo;

            desc = ['EN atypes ',num2str(atype),' ',num2str(atype2),' pca '];
            m1 = Mixer(parIn, mixType, desc, functype);
            m1.fixed = fixed;
            m1.hybrid = 1;
            minfo.mixer = m1;
            minfo.type = 'ENbond';
            mixInfo{end+1} = minfo;
            
            desc = ['E2 atypes ',num2str(atype),' ',num2str(atype2),' pca '];
            m1 = Mixer(parIn, mixType, desc, functype);
            m1.fixed = fixed;
            minfo.mixer = m1;
            minfo.type = 'E2bond';
            mixInfo{end+1} = minfo;

%             parIn = [0 0 0];
%             mixType = 32; % MM stretch mixer
%             desc = ['KE MM atypes ',num2str(atype),' ',num2str(atype2),' pca '];
%             functype = 0; % scale without const
%             m1 = Mixer(parIn, mixType, desc, functype);
%             m1.fixed = [0 0 0];
%             minfo.mixer = m1;
%             minfo.type = 'KEMM';
%             minfo.atype1 = atype;
%             minfo.atype2 = atype2;
%             mixInfo{end+1} = minfo;
            
         end
      end
   end
end

minfo.mixer = Mixer([1 0],4,'e2.HH',2);
minfo.mixer.fixed = [0 1];
minfo.type = 'E2bond';
minfo.atype1 = 1;
minfo.atype2 = 1;
mixInfo{end+1} = minfo;

%% add mixers to all models
allModels = {mtrain{:},mtest{:}};
for imix = 1:length(mixInfo)
   minfo = mixInfo{imix};
   atype1 = minfo.atype1;
   z1 = Context.atypeToZtype(atype1);
   if (z1 == 1)
      hasSP1 = 0;
      types1 = 1;
   else
      hasSP1 = 1;
      types1 = [1 2];
   end
   if (isfield(minfo,'atype2'))
      atype2 = minfo.atype2;
      z2 = Context.atypeToZtype(atype2);
      if (z2 == 1)
         hasSP2 = 0;
         types2 = [1];
      else
         hasSP2 = 1;
         types2 = [1 2];
      end
   end
   for imod = 1:length(allModels)
      mod = allModels{imod};
      switch minfo.type
         case 'KEdiag'
            for itype = types1
               mod.addKEmodDiag(atype1,itype,minfo.mixer);
            end
         case 'KEdiags'
            for itype = types1
               mod.addKEmodDiag(atype1,1,minfo.mixer);
            end
         case 'KEdiagp'
            for itype = types1
               mod.addKEmodDiag(atype1,2,minfo.mixer);
            end
         case 'KE1s1s'
            diag = 1;
            mod.addKEcore(diag,atype1,minfo.mixer);
         case 'KE1s2s'
            diag = 0;
            mod.addKEcore(diag,atype1,minfo.mixer);
         case 'ENdiag'
            for itype = types1
               mod.addENmodDiag(atype1,itype,minfo.mixer);
            end
         case 'E2diag'
               mod.addH2modDiag(atype1,minfo.mixer);
         case 'KEbond'
            mod.addKEmodBondedh(atype1,atype2,minfo.mixer);
         case 'ENbond'
            mod.addENmodBonded1h(atype1,atype2,minfo.mixer);
            if (atype1 ~= atype2)
               mod.addENmodBonded1h(atype2,atype1,minfo.mixer);
            end
%          case 'KEMM'
%             mod.addMMStretch('KE',atype1,atype2, minfo.mixer);
         case 'E2bond'
            mod.addH2modOffDiag(atype1,atype2,minfo.mixer);
         otherwise
            error('Context.makefitme unrecognized mixer type');
      end
   end
end


%% make the fitme classes

ftrain = Fitme;
for imod = 1:length(mtrain);
   if (mtrain{imod}.natom == 4)
      plotnumber = 800;
   else
      plotnumber = 801;
   end
   ftrain.addFrag(mtrain{imod},HLTrain{imod},plotnumber);
end
ftrain.includeKE = 1;
ftrain.includeEN = ones(1,6);
ftrain.includeE2 = 1;
ftrain.silent = 1;
ftrain.setEnvs(envsTrain);
% setEnvs calculates the HL values of everything we are fitting to, so the
% HLs are no longer needed. By removing these from fitme, we make the fitme
% object quite a bit smaller.
ftrain.HLs = [];
ftrain.plot = 1;

ftest = Fitme;
for imod = 1:length(mtest);
   if (mtest{imod}.natom == 4)
      plotnumber = 900;
   else
      plotnumber = 901;
   end
   ftest.addFrag(mtest{imod},HLTest{imod},plotnumber);
end
ftest.includeKE = 1;
ftest.includeEN = ones(1,6);
ftest.includeE2 = 1;
ftest.silent = 1;
ftest.setEnvs(envsTest);
% setEnvs calculates the HL values of everything we are fitting to, so the
% HLs are no longer needed. By removing these from fitme, we make the fitme
% object quite a bit smaller.
ftest.HLs = [];
ftest.plot = 1;
