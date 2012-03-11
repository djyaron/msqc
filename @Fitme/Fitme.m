classdef Fitme < handle
   properties
      models  % {1,nmodels}  cell array of models
      HLs     % {1,nmodels}  cell array to fit to
      HLKE    % {1,nmodels}(1,nenv+1) KE energy
      HLEN    % {1,nmodels}(natom,nenv+1) electron-nuclear interaction
      
      includeKE % include kinetic energy in fit
      includeEN % {1,natoms} include individual electron-nuclear operators
      exactDensity % re-evaluate density matrix on every call to err()
      
      parHF   % Last parameters for which HF was solved
   end
   methods
      function res = Fitme
         res.models = cell(0,0);
         res.HLs    = cell(0,0);
         res.HLKE   = cell(0,0);
         res.exactDensity = 0;
         res.includeKE = 1;
         res.includeEN = [];
      end
      function addFrag(obj,model,HL)
         obj.models{1,end+1} = model;
         obj.HLs{1,end+1} = HL;
         ke = zeros(1,HL.nenv+1);
         for ienv = 0:HL.nenv
            ke(1,ienv+1) = sum(sum( HL.partitionE1(ienv, HL.KE) ));
         end
         obj.HLKE{1,end+1} = ke;
         en = zeros(HL.natom, HL.nenv+1);
         for iatom = 1:HL.natom
            for ienv = 0:HL.nenv
               en(iatom,ienv+1) = sum(sum( HL.partitionE1(ienv, ...
                  HL.H1en(:,:,iatom)) ));
            end
         end
         obj.HLEN{1,end+1} = en;
         if (size(obj.includeEN,1) == 0)
            obj.includeEN = ones(1,HL.natom);
         end
      end
      function res = nmodels(obj)
         res = size(obj.models,2);
      end
      function updateDensity(obj,par)
         for i = 1:obj.nmodels
            obj.models{i}.setPar(par);
            obj.models{i}.solveHF();
         end
         obj.parHF = par;
      end
      function res = err(obj,par)
         flip = 0;
         if (size(par,1)>size(par,2))
            par = par';
            flip = 1;
         end
         disp(['Fitme.err called with par = ',num2str(par)]);
         for imod = 1:obj.nmodels
            obj.models{imod}.setPar(par);
         end
         if (size(obj.parHF,1) == 0)
            dpar = 1;
         else
            dpar = max(abs(obj.parHF-par));
         end
         if (obj.exactDensity && (dpar > -1))
            disp(['solving for density matrices']);
            for imod = 1:obj.nmodels
               obj.models{imod}.solveHF();
            end
            obj.parHF = par;
         end
         % Do sum over all orbitals
         sumRange = cell(1,1);
         ic = 0;
         for imod = 1:obj.nmodels
            sumRange{1,1} = 1:obj.models{imod}.nbasis;
            if (obj.includeKE == 1)
               for ienv = 0:obj.models{imod}.nenv
                  ic = ic + 1;
                  res(ic) = obj.models{imod}.partitionE1(ienv , ...
                     obj.models{imod}.KE(ienv), sumRange) - ...
                     obj.HLKE{1,imod}(ienv+1);
               end
            end
            for iatom = 1:obj.models{imod}.natom
               if (obj.includeEN(iatom) == 1)
                  for ienv = 0:obj.models{imod}.nenv
                     ic = ic+1;
                     res(ic) = obj.models{imod}.partitionE1(ienv, ...
                        obj.models{imod}.H1en(iatom,ienv), sumRange) - ...
                        obj.HLEN{1,imod}(iatom,ienv+1);
                  end
               end
            end
         end
         disp(['RMS err = ',num2str(sqrt(res*res')/ic)]);
         % figure(101)
         % plot(res,'r.');
         if (flip == 1)
            res = res';
         end
      end
      function res = errDiffs(obj,par)
         disp(['Fitme.err called with par = ',num2str(par)]);
         for imod = 1:obj.nmodels
            obj.models{imod}.setPar(par);
         end
         if (size(obj.parHF,1) == 0)
            dpar = 1;
         else
            dpar = max(abs(obj.parHF-par));
         end
         if (obj.exactDensity && (dpar > -1))
            disp(['solving for density matrices']);
            for imod = 1:obj.nmodels
               obj.models{imod}.solveHF();
            end
            obj.parHF = par;
         end
         % Do sum over all orbitals
         sumRange = cell(1,1);
         ic = 0;
         for imod = 1:obj.nmodels
            if (obj.includeKE)
               sumRange{1,1} = 1:obj.models{imod}.nbasis;
               LL0 = obj.models{imod}.partitionE1(0 , ...
                  obj.models{imod}.KE(0), sumRange);
               HL0 =  obj.HLKE{1,imod}(1);
               %ic = ic+1;
               %res(ic) = LL0-HL0;
               for ienv = 1:obj.models{imod}.nenv
                  ic = ic + 1;
                  res(ic) = (obj.models{imod}.partitionE1(ienv , ...
                     obj.models{imod}.KE(ienv), sumRange)-LL0) - ...
                     (obj.HLKE{1,imod}(ienv+1)-HL0);
               end
            end
            for iatom = 1:obj.models{imod}.natom
               if (obj.includeEN(iatom))
                  LL0 = obj.models{imod}.partitionE1(0, ...
                     obj.models{imod}.H1en(iatom), sumRange);
                  HL0 =  obj.HLEN{1,imod}(iatom,1);
                  %ic = ic+1;
                  %res(ic) = LL0-HL0;
                  for ienv = 0:obj.models{imod}.nenv
                     ic = ic+1;
                     res(ic) = (obj.models{imod}.partitionE1(ienv, ...
                        obj.models{imod}.H1en(iatom), sumRange)-LL0) - ...
                        (obj.HLEN{1,imod}(iatom,ienv+1)-HL0);
                  end
               end
            end
         end
         figure(101)
         plot(res,'r.');
         disp(['RMS err = ',num2str(sqrt(res*res')/ic)]);
      end
   end   
end

