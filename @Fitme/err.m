function [res plotnum etype modelnum envnum] = err(obj,par)
flip = 0; % to handle fit routines that pass row or column
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
if (~obj.parallel)
   dpar = obj.updateDensity();
else
   if ((size(obj.parHF,1) == 0) || ...
         (length(obj.parHF) ~= length(par) ) )
      dpar = 1e10;
   else
      dpar = max(abs(obj.parHF-par));
   end
   % need to reset saved densities if changes are bigger than epsDensity2
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
   
   % Write out matlab files that hold each of the models
   maxEnv = 0;
   for imod = 1:obj.nmodels
      fileName = [obj.scratchDir,'fitmeMod',num2str(imod),'.mat'];
      mod = obj.models{imod};
      save(fileName,'mod');
      maxEnv = max([maxEnv, length(obj.envs{1,imod})]);
   end
   
   % create structure saying which calcs to do
   calcs = {};
   calcsInv = zeros(obj.nmodels,maxEnv);
   for imod = 1:obj.nmodels
      envs = obj.envs{1,imod};
      for ienv = 1:length(envs)
         t1.imod = imod;
         t1.ienv = envs(ienv);
         calcs{end+1} = t1;
         calcsInv(imod,ienv) = length(calcs);
      end
   end
   
   % Do the calcs and save results
   ncalc = length(calcs);
   calcRes = cell(ncalc,1);
   % if Etot is fit, need all operators
   if (~obj.includeEtot)
      includeKE = obj.includeKE;
      includeEN = obj.includeEN;
      includeE2 = obj.includeE2;
   else
      includeKE = 1;
      includeEN = ones(1,20);
      includeE2 = 1;
   end
   scratchDir = obj.scratchDir;
   % disp('starting parfor loop');
   parfor icalc = 1:ncalc
      imod = calcs{icalc}.imod;
      ienv = calcs{icalc}.ienv;
      [ke, en, e2, newDensity] = Fitme.modelCalcParallel(imod,ienv,redoDensity,...
         includeKE,includeEN,includeE2,scratchDir);
      t1 = [];
      t1.ke = ke;
      t1.en = en;
      t1.e2 = e2;
      t1.density = newDensity;
      calcRes{icalc} = t1;
   end
   
   for icalc = 1:ncalc
      imod = calcs{icalc}.imod;
      ienv = calcs{icalc}.ienv;
      obj.models{imod}.densitySave{ienv+1} = calcRes{icalc}.density;
   end
end
% to allow control-C to stop the job
pause(0.1);

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
   % hold pointer to start of this model in resTot
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
      if (~obj.parallel)
         modpred = obj.models{imod}.EKE(obj.envs{1,imod});
      else
         nenvs = length(obj.envs{1,imod});
         modpred = zeros(1,nenvs);
         for i=1:nenvs
            icalc = calcsInv(imod,i);
            modpred(i) = calcRes{icalc}.ke;
         end
      end
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
      if (calcEN( obj.models{imod}.Z(iatom) ))
         hlevel = obj.HLEN{imod}(iatom,:);
         if (~obj.parallel)
            modpred = obj.models{imod}.Een(iatom,obj.envs{1,imod});
         else
            nenvs = length(obj.envs{1,imod});
            modpred = zeros(1,nenvs);
            for i=1:nenvs
               icalc = calcsInv(imod,i);
               modpred(i) = calcRes{icalc}.en(iatom);
            end
         end
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
      if (~obj.parallel)
         modpred = obj.models{imod}.E2(obj.envs{1,imod});
      else
         nenvs = length(obj.envs{1,imod});
         modpred = zeros(1,nenvs);
         for i=1:nenvs
            icalc = calcsInv(imod,i);
            modpred(i) = calcRes{icalc}.e2;
         end
      end
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
if (~obj.silent)
   disp(['RMS err/ndata = ',num2str(sqrt(res*res')/ndat), ...
      ' kcal/mol err = ',num2str(sqrt(res*res'/ndat)*627.509)]);
end
obj.itcount = obj.itcount + 1;
obj.errTrain(obj.itcount) = norm(res);
if (size(obj.testFitme,1) > 0)
   err1 = obj.testFitme.err(par);
   obj.errTest(obj.itcount) = norm(err1);
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