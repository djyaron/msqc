classdef Fitme < handle
   properties
      models  % {1,nmodels}  cell array of models
      HLs     % {1,nmodels}  cell array to fit to
      HLKE    % {1,nmodels}(1,nenv) KE energy
      HLEN    % {1,nmodels}(natom,nenv) electron-nuclear interaction
      HLE2    % {1,nmodels}(1,nenv) two-elec enerty
      mixers  % {1,nmixer}   cell
      
      envs      % {1,nmodels} list of environments to include in fit
      includeKE % include kinetic energy in fit
      includeEN % {1,Z} include elec-nuc operators for element Z
      includeE2 % include two-elec energy in fit
      
      parHF   % Last parameters for which HF was solved
      epsDensity % re-evaluate density matrix if par change > eps
      epsDensity2 % re-set density matrices if par change > eps
                  % this is because odd parameters can lead to odd
                  % densities, that should be reset.
      
      plot       % Plot results on every call to err()
      LLKE       % {1,nmodels}(1,nenv) used only for plots
      LLEN       % {1,nmodels}(natom,nenv) used only for plots
      LLE2       % {1,nmodels}(1,nenv) used for plots
      plotNumber % (1,nmodels): number for plot of this model
      plotNumErr % plot number for the error plots (default = 799)
      errCalls   % number of calls to the err function
      testFitme  % fitme object that has the test data
      
      itcount    % counter for the err arrays
      errTrain   % error in train as a function of iteration
      errTest    % error in test set
      
      arms       % (npar,narms): for bandit algorithm
      parallel   % true to run err in parallel
      restartFile % place to save intermediate results
      silent     % suppress all displayed output
      
      hftime
   end
   methods (Static)
      [ke, en, e2, newDensity] = ...
            modelCalcParallel(imod,ienv,updateDensity,...
            includeKE,includeEN,includeE2);
   end
   methods
      function res = Fitme
         res.models = cell(0,0);
         res.HLs    = cell(0,0);
         res.HLKE   = cell(0,0);
         res.HLEN   = cell(0,0);
         res.epsDensity = 0.0;
         res.epsDensity2 = 0.01;
         res.includeKE = 1;
         res.includeEN = zeros(1,6);
         res.includeE2 = 0;
         res.parHF = [];
         res.plot = 1;
         res.plotNumber = [];
         res.plotNumErr = 799;
         res.LLKE = cell(0,0);
         res.LLEN = cell(0,0);
         res.errCalls = 0;
         res.itcount = 0;
         res.parallel = 0;
         res.silent = 0;
      end
      function addMixer(obj, mix)
         add = 1;
         for i=1:size(obj.mixers,2)
            if (mix == obj.mixers{1,i})
               add = 0;
               break;
            end
         end
         if (add == 1)
            obj.mixers{1,end+1} = mix;
         end
      end
      function addMixers(obj,mixersIn)
         % mixes is a cell array (1,:) of mixers
         for i=1:size(mixersIn,2)
            obj.addMixer(mixersIn{i});
         end
      end
      function addFrag(obj,model,HL,plotnumber)
         if (nargin < 4)
            plotnumber = 800;
         end
         obj.models{1,end+1} = model;
         obj.addMixers(model.mixers);
         obj.HLs{1,end+1} = HL;
         obj.plotNumber(1,end+1) = plotnumber;
      end
      function setEnvs(obj,envsIn)
         % currently assumes same environments for every frag/model
         if (~iscell(envsIn))
            obj.envs = cell(1,obj.nmodels);
            for i=1:obj.nmodels
               obj.envs{1,i} = envsIn;
            end
         else
            obj.envs = envsIn;
         end
         obj.HLKE = cell(0,0);
         obj.HLEN = cell(0,0);
         obj.HLE2 = cell(0,0);
         obj.LLKE = cell(0,0);
         obj.LLEN = cell(0,0);
         obj.LLE2 = cell(0,0);
         for imod = 1:obj.nmodels
            envs1 = obj.envs{1,imod};
            HL = obj.HLs{imod};
            % will plot against the STO-3G result
            LL = obj.models{imod}.frag;
            obj.HLKE{1,end+1} = HL.EKE(envs1);
            obj.LLKE{1,end+1} = LL.EKE(envs1);
            nenv = size(envs1,2);
            en = zeros(HL.natom, nenv);
            enl = zeros(LL.natom, nenv);
            for iatom = 1:HL.natom
               en(iatom,:) = HL.Een(iatom,envs1);
               enl(iatom,:) = LL.Een(iatom,envs1);
            end
            obj.HLEN{1,end+1} = en;
            obj.LLEN{1,end+1} = enl;
            obj.HLE2{1,end+1} = HL.E2(envs1);
            obj.LLE2{1,end+1} = LL.E2(envs1);
         end
         obj.parHF = [];
      end
      function res = nmodels(obj)
         res = size(obj.models,2);
      end
      function res = npar(obj)
         res = 0;
         for i=1:size(obj.mixers,2)
            res = res + obj.mixers{i}.npar;
         end
      end
      function res = getPars(obj)
         res = zeros(1,obj.npar);
         ic = 1;
         for i = 1:size(obj.mixers,2)
            mtemp = obj.mixers{1,i};
            np = mtemp.npar;
            if (np > 0)
               res(ic:(ic+np-1)) = mtemp.getPars;
            end
            ic = ic + np;
         end
      end
      function setPars(obj,par)
         % sets parameters, and updates densities
         if (size(par,2) ~= obj.npar)
            error(['Fitme.SetPars called with ',num2str(size(par,2)), ...
               ' parameters when ',num2str(obj.npar),' are needed ']);
         end
         ic = 1;
         for i = 1:size(obj.mixers,2)
            mtemp = obj.mixers{1,i};
            np = mtemp.npar;
            if (np > 0)
               mtemp.setPars( par(ic:(ic+np-1)));
            end
            ic = ic + np;
         end
      end
      function dpar = updateDensity(obj)
         par = obj.getPars;
         if ((size(obj.parHF,1) == 0) || ...
            (length(obj.parHF) ~= length(par) ) )
            dpar = 1e10;
         else
            dpar = max(abs(obj.parHF-par));
         end
         if (dpar > obj.epsDensity2)
            for imod = 1:obj.nmodels
               obj.models{imod}.densitySave = ...
                  cell(1,obj.models{imod}.nenv+1);
            end
         end
         if (dpar > obj.epsDensity)
            if (~obj.parallel)
               if (~obj.silent)
                  disp(['solving for density matrices']);
               end
               for imod = 1:obj.nmodels
                  obj.models{imod}.solveHF(obj.envs{1,imod});
               end
            else
               if (~obj.silent)
                  disp(['parallel solving for density matrices']);
               end
               obj.solveHFparallel;
            end
            obj.parHF = par;
         end
      end
      function res = ndata(obj)
         ic = 0;
         for imod = 1:obj.nmodels
            if (obj.includeKE == 1)
               ic = ic + size(obj.HLKE{1,imod},2);
            end
            for iatom = 1:obj.models{imod}.natom
               if (obj.includeEN( obj.models{imod}.Z(iatom) ))
                  ic = ic + size(obj.HLEN{imod}(iatom,:),2);
               end
            end
            if (obj.includeE2 == 1)
               ic = ic + size(obj.HLE2{1,imod},2);
            end
         end
         res = ic;
      end
      function generateArms(obj,narms,plow,phigh)
         rr = rand(obj.npar,narms);
         obj.arms = plow + (phigh-plow).*rr;
      end
      function res = pullArm(obj,iarm)
         imod = randi(obj.nmodels);
         ienv = randi(length(obj.envs{1,imod}));
         obj.models{1,imod}.setPars( obj.arms(:,iarm) );
         obj.models{1,imod}.solveHF(ienv);
         res = (obj.models{1,imod}.EKE(ienv) - obj.HLKE{1,imod}(ienv)).^2;
         for iatom = 1:obj.models{imod}.natom
            res = res + ...
               (obj.HLEN{imod}(iatom,ienv) - ...
               obj.models{imod}.Een(iatom,ienv)).^2;
         end
         res = sqrt(res);
      end
      function res = randMolError(obj,par)
         % selects a random molecule and environment
         % calculates the error for the parameters in par
         % returns a vector of errors to be minimized
         imod = randi(obj.nmodels);
         ienv = randi(length(obj.envs{1,imod}));
         obj.models{1,imod}.setPars( par );
         obj.models{1,imod}.solveHF(ienv);
         res = obj.models{1,imod}.EKE(ienv) - obj.HLKE{1,imod}(ienv);
         for iatom = 1:obj.models{imod}.natom
            temp =  ...
               obj.HLEN{imod}(iatom,ienv) - ...
               obj.models{imod}.Een(iatom,ienv);
            res = [res , temp];
         end
         res = -1.0 * norm(res);
      end
      function res = armError(obj,iarm)
         res = 0;
         for imod = 1:obj.nmodels
            obj.models{1,imod}.setPars( obj.arms(:,iarm) );
            envir =obj.envs{1,imod};
            obj.models{1,imod}.solveHF(envir);
            res = res + norm(obj.models{1,imod}.EKE(envir) - obj.HLKE{1,imod}).^2;
            for iatom = 1:obj.models{imod}.natom
               res = res + ...
                  norm(obj.HLEN{imod}(iatom,:) - ...
                  obj.models{imod}.Een(iatom,envir)).^2;
            end
         end
         res = sqrt(res);
      end
      function res = normErr(obj,par)
         err = obj.err(par);
         res = norm(err);
      end
   end
end