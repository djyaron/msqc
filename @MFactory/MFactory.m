classdef MFactory < handle
   properties
      policy      % cell array of policies
      atomTypes   % atom types used to construct mixInfo
      mixInfo     % cell array of structures specifying mixers
      mixer       % array of unique mixers
      contextName % cell array {mixer,model} of names of context variables
      context     % cell array {mixer,model,ienv} of context values
   end
   methods (Static)
   end
   methods
      function obj = MFactory
         obj.policy  = cell(0,0);
         obj.mixInfo = cell(0,0);
      end
      function res = addPolicy(obj,varargin)
         validParameters = {{'oper','o'},{'func','f'},{'iatom','i'},...
            {'jatom','j'},'sp',{'context','contexts','c'}};
         t1 = validateInput(varargin,validParameters);
         t2.oper = validatestring(t1.oper,{'KE','EN','E2'});
         t2.func = validatestring(t1.func,{'const','scale','interp'});
         t2.sp   = validatestring(t1.sp, ...
            {'sonly','ponly','hybrid','combine','separate'});
         t2.iatom = t1.iatom;
         if (isfield(t1,'jatom'))
            t2.jatom = t1.jatom;
         end
%          % parse out the contexts
%          str = t1.contexts;
%          cs = {};
%          while (~isempty(str))
%             [cs{end+1},str] = strtok(str);
%          end
         if (isfield(t1,'context'))
            t2.context = t1.context;
         end
         obj.policy{end+1} = t2;
         res = t2;
      end
      function addMixer(obj, mix)
         add = 1;
         for i=1:length(obj.mixer)
            if (mix == obj.mixer{1,i})
               add = 0;
               break;
            end
         end
         if (add == 1)
            obj.mixer{1,end+1} = mix;
         end
      end
   end
end
