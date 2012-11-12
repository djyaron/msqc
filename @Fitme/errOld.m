function [res plotnum etype modelnum envnum] = err(obj,par)
flip = 0; % to handle fit routines that pass row or column
if (size(par,1)>size(par,2))
   par = par';
   flip = 1;
end
if (~obj.silent)
   disp(['Fitme.err called with par = ',num2str(par)]);
end
obj.setPars(par);
dpar = obj.updateDensity();

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
for imod = 1:obj.nmodels
   if (obj.includeKE == 1)
      hlevel = obj.HLKE{1,imod};
      modpred = obj.models{imod}.EKE(obj.envs{1,imod});
      t1 = hlevel - modpred;
      n = size(t1,2);
      res(1,ic:(ic+n-1))= t1;
      plotnum(1,ic:(ic+n-1))= obj.plotNumber(imod);
      modelnum(1,ic:(ic+n-1)) = imod;
      envnum(1,ic:(ic+n-1)) = obj.envs{1,imod};
      etype(1,ic:(ic+n-1))= 1;
      ic = ic + n;
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
      if (obj.includeEN( obj.models{imod}.Z(iatom) ))
         hlevel = obj.HLEN{imod}(iatom,:);
         modpred = obj.models{imod}.Een(iatom,obj.envs{1,imod});
         t1 = hlevel - modpred;
         n = size(t1,2);
         res(1,ic:(ic+n-1)) = t1;
         plotnum(1,ic:(ic+n-1))= obj.plotNumber(imod);
         modelnum(1,ic:(ic+n-1)) = imod;
         envnum(1,ic:(ic+n-1)) = obj.envs{1,imod};
         etype(1,ic:(ic+n-1))= 10 + obj.models{imod}.Z(iatom);
         ic = ic + n;
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
   if (obj.includeE2)
      hlevel = obj.HLE2{1,imod};
      modpred = obj.models{imod}.E2(obj.envs{1,imod});
      t1 = hlevel - modpred;
      n = size(t1,2);
      res(1,ic:(ic+n-1))= t1;
      plotnum(1,ic:(ic+n-1))= obj.plotNumber(imod);
      modelnum(1,ic:(ic+n-1)) = imod;
      envnum(1,ic:(ic+n-1)) = obj.envs{1,imod};
      etype(1,ic:(ic+n-1))= 2;
      ic = ic + n;
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
