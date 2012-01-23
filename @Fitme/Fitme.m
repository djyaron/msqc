classdef Fitme < handle
   properties
      models  % {1,nmodels}  cell array of models
      HLs     % {1,nmodels}  cell array to fit to
      HLKE    % {1,nmodels}(1,nenv+1) KE energy
      HLEN    % {1,nmodels}(natom,nenv+1) electron-nuclear interaction
   end
   methods
      function res = Fitme
         res.models = cell(0,0);
         res.HLs    = cell(0,0);
         res.HLKE   = cell(0,0);
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
      end
      function res = nmodels(obj)
         res = size(obj.models,2);
      end
      function updateDensity(obj,par)
         for i = 1:obj.nmodels
            obj.models{i}.par = par;
            obj.models{i}.solveHF();
         end
      end
      function res = err(obj,par)
         ic = 0;
         % Do sum over all orbitals
         sumRange = cell(1,1);
         for imod = 1:obj.nmodels
            obj.models{imod}.par = par;
            sumRange{1,1} = 1:obj.models{imod}.nbasis;
            for ienv = 0:obj.models{imod}.nenv
               ic = ic + 1;
               res(ic) = obj.models{imod}.partitionE1(ienv , ...
                  obj.models{imod}.KE, sumRange) - ...
                  obj.HLKE{1,imod}(ienv+1);
            end
            for ienv = 0:obj.models{imod}.nenv
               for iatom = 1:obj.models{imod}.natom
                  ic = ic+1;
                  res(ic) = obj.models{imod}.partitionE1(ienv, ...
                     obj.models{imod}.H1en(iatom), sumRange) - ...
                     obj.HLEN{1,imod}(iatom,ienv+1);
               end
            end
         end
      end
   end
   
end

