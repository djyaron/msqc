classdef Mixer < handle
   
   properties
      mixType     % 0 for sigmoidal, 1 for linear, 2 for charge dependent
      par         % (1,npar) current parameters
      desc        % string description
      fixed       % (1,npar) 0 if parameter should be fit, 1 if fixed
      hybrid      % 0 = no hybridization, 1 = sigma mods, 2 = pi mods
      funcType    % 1 = interp 2 = scale 3 = scale with const 
                  % 4 = iterp with const
      index       % field that can be used to store an external unique ID
                  % index is not used internally in this class
   end
   
   methods
      function obj = Mixer(parIn,mixType,desc,funcType)
         if (nargin < 1)
            parIn = 0;
         end
         if (nargin < 2)
            mixType = 1;
         end
         if (nargin < 3)
            desc = ' ';
         end
         if (nargin < 4)
            funcType = 1;
         end
         obj.par = parIn;
         obj.mixType = mixType;
         obj.fixed = zeros(size(parIn));
         obj.desc = desc;
         obj.funcType = funcType;
         obj.hybrid = 0;
      end
      function res = deepCopy(obj)
         res = Mixer(obj.par,obj.mixType,obj.desc,obj.funcType);
         res.fixed = obj.fixed;
         res.hybrid = obj.hybrid;
         res.index = obj.index;
      end
      function res = constructionData(obj)
         % data needed to reconstruct this mixer
         res.mixType = obj.mixType;
         res.par = obj.par;
         res.desc = obj.desc;
         res.fixed = obj.fixed;
         res.hybrid = obj.hybrid;
         res.funcType = obj.funcType;
      end
   end
   methods (Static)
      function res = createFromData(dat)
         res = Mixer(dat.par,dat.mixType,dat.desc,dat.funcType);
         res.hybrid = dat.hybrid;
         res.fixed = dat.fixed;
      end
   end
   methods
      function res = npar(obj)
         res = sum(obj.fixed == 0);
      end
      function setPars(obj, pars)
         ic = find(obj.fixed == 0);
         if (size(ic,2) > 0)
            obj.par(ic) = pars;
         end
      end
      function res = getPars(obj)
         ic = find(obj.fixed == 0);
         if (~isempty(ic))
            res = obj.par(ic);
         else
            res = [];
         end
      end
      function res = print(obj)
         types = {'sigmoid','linear','ch-dep','bo-dep','bl-dep','bo-bl-dep'};
         types{22 + 1} = 'ch-dep-quad';
         types{11 + 1} = 'context-atom';
         types{12 + 1} = 'context-bond';
         types{32 + 1} = 'MM stretch';
         ftypes = {' ','mult','mult-c','mix-c'};
         res = [obj.desc,' ',types{obj.mixType+1},' ',...
            ftypes{obj.funcType}];
         if (obj.hybrid == 1)
            res = [res,' sig'];
         end
         if (obj.hybrid == 2)
            res = [res,' pi'];
         end
         for i = 1:size(obj.par,2)
            res = [res,' ',num2str(obj.par(i))];
            if (obj.fixed(i))
               res = [res,'*'];
            end
         end
         disp(res);
      end
   end
end

