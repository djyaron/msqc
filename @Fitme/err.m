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

nmodels = obj.nmodels;
ncalc = length([obj.envs{:}]);

% Decide how to distribute the work among the pool of workers.
minWorkload = 10;  % Want many jobs per worker.
if (obj.parallel == 0)
   jobsPerGroup = ncalc;
else
   jobsPerGroup = max([1 ceil(ncalc/(minWorkload * obj.parallel))]);
end
groupSize = ceil(ncalc/jobsPerGroup);

calcRes = cell(jobsPerGroup, groupSize);
calcs = cell(jobsPerGroup, groupSize);

% Write out matlab files that hold each of the models.
maxEnv = 0;
if (~isempty(obj.scratchDir) && exist(obj.scratchDir, 'dir') ~= 7)
   mkdir(obj.scratchDir);
end
for imod = 1:nmodels
   fileName = [obj.scratchDir,'fitmeMod',num2str(imod),'.mat'];
   mod = obj.models{imod}; %#ok<NASGU>
   save(fileName,'mod');
   maxEnv = max([maxEnv, length(obj.envs{1,imod})]);
end

% Create structure saying which calcs to do.
calcsInv = zeros(nmodels, maxEnv);
icalc = 1;
for imod = 1:obj.nmodels
   envs = obj.envs{1,imod};
   for ienv = 1:length(envs)
      t1.imod = imod;
      t1.ienv = envs(ienv);
      calcs{icalc} = t1;
      calcsInv(imod, ienv) = icalc;
      icalc = icalc + 1;
   end
end

% Do the calcs and save results.
% If Etot is fit, need all operators.
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
parfor j = 1:groupSize
   lastmod = -1;
   for i = 1:jobsPerGroup
      if (~isempty(calcs{i, j}))
         imod = calcs{i, j}.imod;
         ienv = calcs{i, j}.ienv;
         if (imod ~= lastmod)
            % Retrieve model from save file.
            fileName = [scratchDir,'fitmeMod',num2str(imod),'.mat'];
            saveData = load(fileName);
            data = saveData.mod;
         end
         lastmod = imod;
         [ke, en, e2, newDensity, orb, Eorb, Ehf] = ...
            data.modelCalcParallel(ienv,redoDensity,...
            includeKE,includeEN,includeE2);
         t1 = struct;
         t1.ke = ke;
         t1.en = en;
         t1.e2 = e2;
         t1.density = newDensity;
         t1.orb  = orb;
         t1.Eorb = Eorb;
         t1.Ehf  = Ehf;
         calcRes{i, j} = t1;
      end
   end
end

if (redoDensity)
   for icalc = 1:ncalc
      imod = calcs{icalc}.imod;
      ienv = calcs{icalc}.ienv;
      obj.models{imod}.densitySave{ienv+1} = calcRes{icalc}.density;
      if (ienv == 0)
         obj.models{imod}.orb  = calcRes{icalc}.orb;
         obj.models{imod}.Eorb = calcRes{icalc}.Eorb;
         obj.models{imod}.Ehf  = calcRes{icalc}.Ehf;
      else
         obj.models{imod}.orbEnv(:,:,ienv) = calcRes{icalc}.orb;
         obj.models{imod}.EorbEnv(:,ienv)  = calcRes{icalc}.Eorb;
         obj.models{imod}.EhfEnv(1,ienv)   = calcRes{icalc}.Ehf;
      end
   end
end

% To allow control-C to stop the job.
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
      nenvs = length(obj.envs{1,imod});
      modpred = zeros(1,nenvs);
      for i=1:nenvs
         icalc = calcsInv(imod,i);
         modpred(i) = calcRes{icalc}.ke;
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
      if (calcEN(obj.models{imod}.Z(iatom)))
         hlevel = obj.HLEN{imod}(iatom,:);
         nenvs = length(obj.envs{1,imod});
         modpred = zeros(1,nenvs);
         for i=1:nenvs
            icalc = calcsInv(imod,i);
            modpred(i) = calcRes{icalc}.en(iatom);
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
      nenvs = length(obj.envs{1,imod});
      modpred = zeros(1,nenvs);
      for i=1:nenvs
         icalc = calcsInv(imod,i);
         modpred(i) = calcRes{icalc}.e2;
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