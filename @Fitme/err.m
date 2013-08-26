function [res, plotnum, etype, modelnum, envnum] = err(obj, par)
% Description:
%   Gateway error routine.
%
% Input:
%   par:
%
% Output:
%   res:
%   plotnum:
%   etype: Operator type. 1->KE, 2->E2, 3->Etot, 10+#->EN for #.
%   modelnum:
%   envnum:

if (nargin < 2)
   par = obj.getPars;
end

flip = 0;  % To handle fit routines that pass row or column.
if (size(par,1)>size(par,2))
   par = par';
   flip = 1;
end
if (~obj.silent)
   disp(['Fitme.err called with par = ',num2str(par)]);
end
if (~isempty(obj.cset))
   obj.cset.saveIndices;
end

obj.setPars(par);
if ((size(obj.parHF,1) == 0) || ...
      (length(obj.parHF) ~= length(par) ) )
   dpar = 1e10;
else
   dpar = max(abs(obj.parHF-par));
end

% Need to reset saved densities if changes are bigger than epsDensity2.
if (dpar > obj.epsDensity2)
   message1 = 'Densities reset';
   for imod = 1:obj.nmodels
      obj.models{imod}.densitySave = ...
         cell(1,obj.models{imod}.nenv+1);
   end
else
   message1 = '';
end

if (dpar > obj.epsDensity)
   if (~obj.silent)
      disp(['density matrices will be recalculated ',message1]);
   end
   redoDensity = 1;
else
   redoDensity = 0;
end

% How many H1,H2 matrices to save before entering the parfor loop
obj.parallel = matlabpool('size');
if (obj.parallel == 0)
   for imod = 1:obj.nmodels
      envs = obj.envs{1,imod};
      for ienv = 1:length(envs)
         obj.models{imod}.solveHF(obj.envs{1,imod}, ...
            obj.HFeps,obj.HFmaxit,obj.HFminit);
      end
   end
else
   saveForParallel = 20 * obj.parallel;
   whichCalc = cell(1,saveForParallel);
   H1p  = cell(1,saveForParallel);
   H2p  = cell(1,saveForParallel);
   Sp  = cell(1,saveForParallel);
   Enucp  = cell(1,saveForParallel);
   Nelecp  = cell(1,saveForParallel);
   guessDensity = cell(1,saveForParallel);
   orbp  = cell(1,saveForParallel);
   Eorbp = cell(1,saveForParallel);
   Ehfp  = cell(1,saveForParallel);
   HFeps = obj.HFeps;
   HFmaxit = obj.HFmaxit;
   HFminit = obj.HFminit;
   istore = 0;
   for imod = 1:obj.nmodels
      mod1 = obj.models{imod};
      envs = obj.envs{1,imod};
      for indexEnv = 1:length(envs)
         ienv = envs(indexEnv);
         istore = istore +1;
         whichCalc{istore}.imod = imod;
         whichCalc{istore}.ienv = ienv;
         H1p{istore} = obj.models{imod}.H1(ienv);
         H2p{istore} = obj.models{imod}.H2(ienv);
         Sp{istore} = obj.models{imod}.S;
         Enucp{istore} = obj.models{imod}.Hnuc(ienv);
         Nelecp{istore} = obj.models{imod}.frag.nelec;
         % copied directly from Model3.HartreeFock
         if ((size(mod1.densitySave{ienv+1},1) == 0) && ...
               (size(mod1.densitySave{1},1) == 0))
            guessDensity{istore} = mod1.frag.density(ienv);
         elseif (size(mod1.densitySave{ienv+1},1) == 0)
            guessDensity{istore} = mod1.densitySave{1};
         else
            guessDensity{istore} = mod1.densitySave{ienv+1};
         end
         % do calculations if we reach saveForParallel, loop is terminating
         if (istore == saveForParallel)
            doCalcs = saveForParallel;
         elseif ((imod == obj.nmodels) && (indexEnv == length(envs)))
            doCalcs = istore;
         else
            doCalcs = 0;
         end
         if (doCalcs > 0)
            parfor ip = 1:doCalcs
               [orbp{ip},Eorbp{ip},Ehfp{ip}] = parallelHF(H1p{ip},H2p{ip} ...
                  ,Sp{ip},Enucp{ip},Nelecp{ip} ...
                  ,guessDensity{ip},HFeps,HFmaxit,HFminit);
            end
            for ip = 1:doCalcs
               smod = obj.models{whichCalc{ip}.imod};
               senv = whichCalc{ip}.ienv;
               if (senv == 0)
                  smod.orb = orbp{ip};
                  smod.Eorb = Eorbp{ip};
                  smod.Ehf = Ehfp{ip};
               else
                  smod.orbEnv(:,:,senv) = orbp{ip};
                  smod.EorbEnv(:,senv) = Eorbp{ip};
                  smod.EhfEnv(1,senv) = Ehfp{ip};
               end
               smod.densitySave{senv+1} = smod.density(senv);
            end
            istore = 0;
         end
      end
   end
end
% To allow control-C to stop the job.
pause(0.01);

doPlots = obj.plot && (dpar > 1.0e-4);

if (doPlots)
   for i=unique([obj.plotNumber])
      figure(i);
      clf;
   end
end

ic = 1;
ndat = obj.ndata;
res = zeros(1,ndat);
plotnum = zeros(1,ndat);
etype = zeros(1,ndat);
if (obj.includeEtot)
   calcKE = 1;
   calcEN = ones(1,20);
   calcE2 = 1;
   % Hold pointer to start of this model in resTot.
   rangeMod = cell(1,obj.nmodels);
   ntot = 0;
   for imod = 1:obj.nmodels
      rangeMod{imod} = (ntot+1):(ntot+length(obj.envs{1,imod}));
      ntot = ntot + length(obj.envs{1,imod});
   end
   resTot = zeros(1,ntot);
else
   calcKE = obj.includeKE;
   calcEN = obj.includeEN;
   calcE2 = obj.includeE2;
end
for imod = 1:obj.nmodels
   if (calcKE)
      hlevel = obj.HLKE{1,imod};
      modpred = obj.models{imod}.EKE(obj.envs{imod});
      t1 = hlevel - modpred;
      n = size(t1,2);
      if (obj.includeKE)
         res(1,ic:(ic+n-1))= t1;
         plotnum(1,ic:(ic+n-1))= obj.plotNumber(imod);
         modelnum(1,ic:(ic+n-1)) = imod;
         envnum(1,ic:(ic+n-1)) = obj.envs{1,imod};
         etype(1,ic:(ic+n-1))= 1;
         ic = ic + n;
      end
      if (obj.includeEtot)
         resTot(rangeMod{imod}) = t1;
      end
      if (doPlots)
         figure(obj.plotNumber(imod));
         subplot(4,2,1);
         hold on;
         llevel = obj.LLKE{1,imod};
         plot(llevel,llevel,'k.');
         plot(llevel,hlevel,'r.');
         plot(llevel,modpred,'b.');
         %title('Kinetic E: LL(black) HL(red) model(blue)');
         %xlabel('LL')
         subplot(4,2,2);
         hold on;
         x1 = min(hlevel);
         x2 = max(hlevel);
         plot(hlevel,modpred,'g.');
         plot([x1 x2],[x1 x2],'k-');
         %title('Kinetic E: HL(black) model(red)');
         %xlabel('HL')
      end
   end
   for iatom = 1:obj.models{imod}.natom
      if (calcEN(obj.models{imod}.Z(iatom)))
         hlevel = obj.HLEN{imod}(iatom,:);
         modpred = obj.models{imod}.Een(iatom,obj.envs{imod});
         t1 = hlevel - modpred;
         if (obj.includeEN(obj.models{imod}.Z(iatom) ))
            n = size(t1,2);
            res(1,ic:(ic+n-1)) = t1;
            plotnum(1,ic:(ic+n-1))= obj.plotNumber(imod);
            modelnum(1,ic:(ic+n-1)) = imod;
            envnum(1,ic:(ic+n-1)) = obj.envs{1,imod};
            etype(1,ic:(ic+n-1))= 10 + obj.models{imod}.Z(iatom);
            ic = ic + n;
         end
         if (obj.includeEtot)
            resTot(rangeMod{imod}) = resTot(rangeMod{imod}) + t1;
         end
         if (doPlots)
            if (obj.models{imod}.Z(iatom) == 1)
               frame1 = 3;
               frame2 = 4;
               element = 'H';
            else
               frame1 = 5;
               frame2 = 6;
               element = 'C';
            end
            subplot(4,2,frame1);
            hold on;
            llevel = obj.LLEN{1,imod}(iatom,:);
            plot(llevel,llevel,'k.');
            plot(llevel,hlevel,'r.');
            plot(llevel,modpred,'b.');
            %title(['EN for ',element]);
            %xlabel('LL');
            subplot(4,2,frame2);
            hold on;
            x1 = min(hlevel);
            x2 = max(hlevel);
            plot(hlevel,modpred,'g.');
            plot([x1 x2],[x1 x2],'k-');
            %title(['EN for ',element]);
            %xlabel('HL');
         end
      end
   end
   if (calcE2)
      hlevel = obj.HLE2{1,imod};
      modpred = obj.models{imod}.E2fast(obj.envs{imod});
      t1 = hlevel - modpred;
      if (obj.includeE2)
         n = size(t1,2);
         res(1,ic:(ic+n-1))= t1;
         plotnum(1,ic:(ic+n-1))= obj.plotNumber(imod);
         modelnum(1,ic:(ic+n-1)) = imod;
         envnum(1,ic:(ic+n-1)) = obj.envs{1,imod};
         etype(1,ic:(ic+n-1))= 2;
         ic = ic + n;
      end
      if (obj.includeEtot)
         resTot(rangeMod{imod}) = resTot(rangeMod{imod}) + t1;
      end
      if (doPlots)
         figure(obj.plotNumber(imod));
         subplot(4,2,7);
         hold on;
         llevel = obj.LLE2{1,imod};
         plot(llevel,llevel,'k.');
         plot(llevel,hlevel,'r.');
         plot(llevel,modpred,'b.');
         %title('E2: LL(black) HL(red) model(blue)');
         %xlabel('LL')
         subplot(4,2,8);
         hold on;
         x1 = min(hlevel);
         x2 = max(hlevel);
         plot(hlevel,modpred,'g.');
         plot([x1 x2],[x1 x2],'k-');
         %title('E2: HL(black) model(red)');
         %xlabel('HL')
      end
   end
end
if (obj.includeEtot)
   for imod = 1:obj.nmodels
      t1 = resTot(rangeMod{imod});
      n = size(t1,2);
      res(1,ic:(ic+n-1))= t1;
      plotnum(1,ic:(ic+n-1))= obj.plotNumber(imod);
      modelnum(1,ic:(ic+n-1)) = imod;
      envnum(1,ic:(ic+n-1)) = obj.envs{1,imod};
      etype(1,ic:(ic+n-1))= 3;
      ic = ic + n;
   end
end
% For backwards compatability with fitme saved before this variable
% was added.
if (isempty(obj.cost))
   obj.cost = 0;
end
if (obj.cost > 0)
   obj.setCostVector;
   res(ic) = (ndat * obj.cost/627.509) * norm(obj.costVector(:) .* par(:));
   plotnum(ic) = 0;
   modelnum(ic) = 0;
   envnum(ic) = -1;
   etype(ic) = 0;
   ic = ic + 1;
end
if (~obj.silent)
   if (obj.cost > 0)
      nc = res(1:(ic-2));
      wc = res(ic-1);
      disp(['RMS err/ndata = ',num2str(sqrt(nc*nc')/ndat), ...
         ' kcal/mol err = ',num2str(sqrt(nc*nc'/ndat)*627.509) ...
         ' kcal/mol cost = ',num2str(wc*627.509/ndat)]);
   else
      disp(['RMS err/ndata = ',num2str(sqrt(res*res')/ndat), ...
         ' kcal/mol err = ',num2str(sqrt(res*res'/ndat)*627.509)]);
   end
end
obj.itcount = obj.itcount + 1;
obj.errTrain(obj.itcount) = norm(res);
if (size(obj.testFitme,1) > 0)
   err1 = obj.testFitme.err(par);
   obj.errTest(obj.itcount) = norm(err1);
end
if (~isempty(obj.operWeights))
   % Has fields KE, EN(1...Zmax), E2 and Etot.
   ike = (etype==1);
   res(ike) = res(ike) * obj.operWeights.KE;
   % Find all atom types.
   atypes = unique(etype(etype>10));
   for atype = atypes
      iz = (etype==atype);
      res(iz) = res(iz) * obj.operWeights.EN(atype-10);
   end
   i2 = (etype==2);
   res(i2) = res(i2)*obj.operWeights.E2;
   itot = (etype==3);
   res(itot) = res(itot) * obj.operWeights.Etot;
   if (~obj.silent)
      disp(['Weighted RMS err/ndata = ',num2str(sqrt(res*res')/ndat), ...
         ' kcal/mol err = ',num2str(sqrt(res*res'/ndat)*627.509)]);
   end
end

if (doPlots)
   figure(obj.plotNumErr);
   if (obj.errCalls == 0)
      hold off;
      title('log10(error) for test (red+) and train (blue o)');
   else
      hold on;
   end
   plot(obj.errCalls+1, log10(norm(res)/length(res)),'bo');
   if (size(obj.testFitme,1) > 0)
      hold on;
      plot(obj.errCalls+1, log10(norm(err1)/length(err1)),'r+');
   end
end
if (flip == 1)
   res = res';
end
obj.errCalls = obj.errCalls + 1;

if (dpar > 1.0e-4)
   if (~isempty(obj.restartFile))
      if (~obj.silent)
         disp('saving restart file');
      end
      ptSave = par;
      itSave = obj.itcount;
      errTrainSave = obj.errTrain;
      errTestSave = obj.errTest;
      save(obj.restartFile,'ptSave','itSave', ...
         'errTrainSave','errTestSave');
   end
end

end