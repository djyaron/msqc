classdef Context < handle
   % intended use to get context from a train and test set
   % mtrain is a set of models to train on, and mtest is the test set
   %    iatom = 1
   %    c1 = Context(length(mtrain))
   %    for i =1:length(mtrain)
   %       c1.addModel(mtrain{i},iatom) % add jatom to list for bond
   %    end
   %    c1.extractFeatures;
   %    % features for any model and atom can be obtained from
   %    c1.project(mod1, iatom);
   properties
      ndim       % length of vectors for feature extraction
      data       % (nmod, ndim)
   end
   
   methods (Static)
      fillInContexts(modsTrain,envsTrain,modsTest,envsTest)
      % modsTrain and modsTest are cell arrays of models
      % This function fills in all the necessary atom and bond contexts
   end % static methods
   methods
      function obj = Context(nmod)
         % preconfigures the storage assuming addModel will be called
         % n mod times
         ndim = 12;
         obj.data = zeros(nmod,ndim);
      end
      function res = analyzeBond(obj,mod,ienv,iatom,jatom)
         % Determine variables that characterize a bond from iatom to jatom
      end
      function addModel(obj,mod,ienv,iatom,jatom)
         % creates a unique ordering for the bond values, for a particular
         % atom (if only iatom passed) or bond (if both iatom and jatom are
         % passed), and adds these to a list for feature extraction
      end
      function extractFeatures(obj)
         % performs feature extraction on data that has been added
      end
      function project(obj,mod,ienv,iatom,jatom)
         % projects an atom or bond onto the features
      end
   end % methods
end

