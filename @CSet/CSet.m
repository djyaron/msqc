classdef CSet < handle
   % Set of contexts corresponding to set of models and mixers
   % The context is a function of both the model and the mixer. The index
   % property of the models and mixers is used as a unique key identifying
   % the object. 
   % For diagonal mixers, the context is a function of the atom number and
   % environment number. This is stored internally as a {natom,nenv}
   % cell array
   % For offdiagonal mixers, teh context is a function of the bonded atoms.
   % This is stored internally as a {natom,natom,nenv} cell array
   properties
      models       % list of models
      mixers       % list of mixer indices
      context      % cell array of contexts {modelIndex, mixerIndex}
                   % for atoms, it is a matrix of size 
                   % (nc,natom,nenv) for atoms
                   % (nc,natom,natom,nenv} for bonds
   end
   methods
      function res = CSet
         res.models = cell(0,0);
         res.mixers = cell(0,0);
      end
      function [res added] = modelIndex(obj,model)
         % return index of model, adding it if needed
         for i = 1:length(obj.models)
            if (model == obj.models{i})
               res = i;
               added = 0;
               return;
            end
         end
         obj.models{end+1} = model;
         added = 1;
         res = length(obj.models);
      end
      function [res added] = mixerIndex(obj,mixer)
         % return index of mixer, adding it if needed
         for i = 1:length(obj.mixers)
            if (mixer == obj.mixers{i})
               res = i;
               added = 0;
               return;
            end
         end
         obj.mixers{end+1} = mixer;
         added = 1;
         res = length(obj.mixers);
      end
      function saveIndices(obj)
         for imod = 1:length(obj.models)
            obj.models{imod}.index = imod;
         end
         for imix = 1:length(obj.mixers)
            obj.mixers{imix}.index = imix;
            obj.mixers{imix}.cset = obj;
         end
      end
      function res = getContext(obj,model,mixer,iatom,jatom,ienv)
         %imod = obj.modelIndex(model);
         %imix = obj.mixerIndex(mixer);
         %if ((imix > size(obj.context,2)) || (imod>size(obj.context,1)))
         %   disp('what?');
         %end
         %t1 = obj.context{imod,imix};
         t1 = obj.context{model.index,mixer.index};
         if (iatom == jatom)
            res = t1(:,iatom,ienv+1);
         else
            res = t1(:,iatom,jatom,ienv+1);
         end
      end
   end
end
