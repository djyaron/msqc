classdef MSet < handle
   properties
      models  % cell array of current models
      envs    % cell array of the envs to include for each model
      pnum    % number used to specify plot, or printEdetails
      HLfrag  % cell array of high level frags corresponding to model
   end
   methods (Static)
      function res = fromFitme(f1)
         res = MSet;
         res.models = f1.models;
         res.envs = f1.envs;
         res.pnum = f1.plotNumber;
      end
   end
   methods
      function obj = MSet
         obj.models  = cell(0,0);
         obj.envs   = cell(0,0);
         obj.HLfrag = cell(0,0);
         obj.pnum = [];
      end
      function addData(obj,fileName,modelsIn,envs,ihl,pnum)
         load(fileName,'LL','HL');
         for i = modelsIn
            obj.models{end+1} = Model3(LL{i,1},LL{i,1},LL{i,1});
            obj.HLfrag{end+1} = HL{i,ihl};
            if (iscell(envs))
               obj.envs{1,end+1} = envs(i);
            else
               obj.envs{1,end+1} = envs;
            end
            obj.pnum(end+1,1) = pnum;
         end
      end
      function res = deepCopy(obj)
         res = MSet;
         for i = 1:length(obj.models)
            m1 = obj.models{i};
            res.models{end+1} = Model3(m1.frag,m1.fnar,m1.fdif);
            res.HLfrag{end+1} = obj.HLfrag{i};
            res.envs{end+1} = obj.envs{i};
            res.pnum(end+1,1) = obj.pnum(i);
         end
      end
      function addSet(obj,newSet)
         for i = 1:length(newSet.models)
            obj.models{end+1} = newSet.models{i};
            obj.HLfrag{end+1} = newSet.HLfrag{i};
            obj.envs{end+1} = newSet.envs{i};
            obj.pnum(end+1,1) = newSet.pnum(i);
         end
      end
      function res = atomTypes(obj)
         allTypes = [];
         for i=1:length(obj.models)
            allTypes = [allTypes,obj.models{i}.aType];
         end
         res = unique(allTypes);
      end
   end
end
