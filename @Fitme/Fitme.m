classdef Fitme < handle
   properties
      models  % {1,nmodels}  cell array of models
      HLs     % {1,nmodels}  cell array to fit to
      HLH1    % {nmodles,1}(1,nenv+1) total H1 energy
   end
   methods
      function res = Fitme
         res.models = cell(0,0);
         res.HLs    = cell(0,0);
         res.HLH1   = cell(0,0);
      end
      function addFrag(obj,model,HL)
         obj.models{1,end+1} = model;
         obj.HLs{1,end+1} = HL;
         h1 = zeros(1,HL.nenv+1);
         for ienv = 0:HL.nenv
            h1(1,ienv+1) = sum(sum( HL.partitionE1(ienv) ));
         end
         obj.HLH1{1,end+1} = h1;
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
         for i = 1:obj.nmodels
            obj.models{i}.par = par;
            for ienv = 0:obj.models{i}.nenv
               ic = ic + 1;
               res(ic) = sum(sum( obj.models{i}.partitionE1(ienv) ))...
                  - obj.HLH1{i}(ienv+1);
            end
         end
      end
   end
   
end

