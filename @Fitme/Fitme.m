classdef Fitme < handle
   properties
      models  % {1,nmodels}  cell array of models
      HLs     % {1,nmodels}  cell array to fit to
      HLKE    % {1,nmodels}(1,nenv) KE energy
      HLEN    % {1,nmodels}(natom,nenv) electron-nuclear interaction
      mixers  % {1,nmixer}   cell
      
      envs      % {1,nmodels} list of environments to include in fit
      includeKE % include kinetic energy in fit
      includeEN % {1,Z} include elec-nuc operators for element Z
      
      parHF   % Last parameters for which HF was solved
      epsDensity % re-evaluate density matrix if par change > eps
      
      plot       % Plot results on every call to err()
      LLKE       % {1,nmodels}(1,nenv) used only for plots
      LLEN       % {1,nmodels}(natom,nenv) used only for plots
      plotNumber % (1,nmodels): number for plot of this model
      plotNumErr % plot number for the error plots (default = 799)
      errCalls   % number of calls to the err function
      testFitme  % fitme object that has the test data
      
      itcount    % counter for the err arrays
      errTrain   % error in train as a function of iteration
      errTest    % error in test set
      
      arms       % (npar,narms): for bandit algorithm
   end
   methods
      function res = Fitme
         res.models = cell(0,0);
         res.HLs    = cell(0,0);
         res.HLKE   = cell(0,0);
         res.HLEN   = cell(0,0);
         res.epsDensity = 0.0;
         res.includeKE = 1;
         res.includeEN = zeros(1,6);
         res.parHF = [];
         res.plot = 1;
         res.plotNumber = [];
         res.plotNumErr = 799;
         res.LLKE = cell(0,0);
         res.LLEN = cell(0,0);
         res.errCalls = 0;
         res.itcount = 0;
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
         obj.envs = cell(1,obj.nmodels);
         for i=1:obj.nmodels
            obj.envs{1,i} = envsIn;
         end
         obj.HLKE = cell(0,0);
         obj.HLEN = cell(0,0);
         obj.LLKE = cell(0,0);
         obj.LLEN = cell(0,0);
         for imod = 1:obj.nmodels
            envs1 = obj.envs{1,i};
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
         if (size(obj.parHF,1) == 0)
            dpar = 1e10;
         else
            dpar = max(abs(obj.parHF-par));
         end
         if (dpar > obj.epsDensity)
            disp(['solving for density matrices']);
            for imod = 1:obj.nmodels
               obj.models{imod}.solveHF(obj.envs{1,imod});
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
            for iatom = 1:obj.HLs{imod}.natom
               if (obj.includeEN( obj.HLs{imod}.Z(iatom) ))
                  ic = ic + size(obj.HLEN{imod}(iatom,:),2);
               end
            end
         end
         res = ic;
      end
      function res = err(obj,par)
         flip = 0; % to handle fit routines that pass row or column
         if (size(par,1)>size(par,2))
            par = par';
            flip = 1;
         end
         disp(['Fitme.err called with par = ',num2str(par)]);
         obj.setPars(par);
         dpar = obj.updateDensity();
         
         doPlots = obj.plot && (dpar > 1.0e-4);
         
         if (doPlots)
            for i=unique(obj.plotNumber)
               figure(i);
               clf;
            end
         end
         
         ic = 1;
         ndat = obj.ndata;
         res = zeros(1,ndat);
         for imod = 1:obj.nmodels
            if (obj.includeKE == 1)
               hlevel = obj.HLKE{1,imod};
               modpred = obj.models{imod}.EKE(obj.envs{1,imod});
               t1 = hlevel - modpred;
               n = size(t1,2);
               res(1,ic:(ic+n-1))= t1;
               ic = ic + n;
               if (doPlots)
                  figure(obj.plotNumber(imod));
                  subplot(3,2,1);
                  hold on;
                  llevel = obj.LLKE{1,imod};
                  plot(llevel,llevel,'k.');
                  plot(llevel,hlevel,'r.');
                  plot(llevel,modpred,'b.');
                  title('Kinetic E: LL(black) HL(red) model(blue)');
                  xlabel('LL')
                  subplot(3,2,2);
                  hold on;
                  x1 = min(hlevel);
                  x2 = max(hlevel);
                  plot(hlevel,modpred,'g.');
                  plot([x1 x2],[x1 x2],'k-');
                  title('Kinetic E: HL(black) model(red)');
                  xlabel('HL')
               end
            end
            for iatom = 1:obj.HLs{imod}.natom
               if (obj.includeEN( obj.HLs{imod}.Z(iatom) ))
                  hlevel = obj.HLEN{imod}(iatom,:);
                  modpred = obj.models{imod}.Een(iatom,obj.envs{1,imod});
                  t1 = hlevel - modpred;
                  n = size(t1,2);
                  res(1,ic:(ic+n-1)) = t1;
                  ic = ic + n;
                  if (doPlots)
                     if (obj.HLs{imod}.Z(iatom) == 1)
                        frame1 = 3;
                        frame2 = 4;
                        element = 'H';
                     else
                        frame1 = 5;
                        frame2 = 6;
                        element = 'C';
                     end
                     subplot(3,2,frame1);
                     hold on;
                     llevel = obj.LLEN{1,imod}(iatom,:);
                     plot(llevel,llevel,'k.');
                     plot(llevel,hlevel,'r.');
                     plot(llevel,modpred,'b.');
                     title(['EN for ',element]);
                     xlabel('LL');
                     subplot(3,2,frame2);
                     hold on;
                     x1 = min(hlevel);
                     x2 = max(hlevel);
                     plot(hlevel,modpred,'g.');
                     plot([x1 x2],[x1 x2],'k-');
                     title(['EN for ',element]);
                     xlabel('HL');
                  end
               end
            end
         end
         disp(['RMS err/ndata = ',num2str(sqrt(res*res')/ndat)]);
         
         if (doPlots)
            figure(obj.plotNumErr);
            if (obj.errCalls == 0)
               hold off;
               title('log10(error) for test (red+) and train (blue o)');
            else
               hold on;
            end
            obj.itcount = obj.itcount + 1;
            plot(obj.errCalls+1, log10(norm(res)/length(res)),'bo');
            obj.errTrain(obj.itcount) = norm(res);
            if (size(obj.testFitme,1) > 0)
               disp('**** TEST SET START ****');
               err1 = obj.testFitme.err(par);
               disp('**** TEST SET END ****');
               hold on;
               plot(obj.errCalls+1, log10(norm(err1)/length(err1)),'r+');
               obj.errTest(obj.itcount) = norm(err1);
            end
         end
         if (flip == 1)
            res = res';
         end
         obj.errCalls = obj.errCalls + 1;
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
         for iatom = 1:obj.HLs{imod}.natom
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
         for iatom = 1:obj.HLs{imod}.natom
            temp =  ...
               obj.HLEN{imod}(iatom,ienv) - ...
               obj.models{imod}.Een(iatom,ienv);
            res = [res , temp];
         end
      end
      function res = armError(obj,iarm)
         res = 0;
         for imod = 1:obj.nmodels
            obj.models{1,imod}.setPars( obj.arms(:,iarm) );
            envir =obj.envs{1,imod};
            obj.models{1,imod}.solveHF(envir);
            res = res + norm(obj.models{1,imod}.EKE(envir) - obj.HLKE{1,imod}).^2;
            for iatom = 1:obj.HLs{imod}.natom
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