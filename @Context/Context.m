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
      data       % (ndim,ndata)
      idata      % how much data is currently filled in
      isBondContext % bool, 1 if the context is for a bond
      mu
      sigma
      coeff
      score
      latent
   end
   
   methods (Static)
      % modsTrain and modsTest are cell arrays of models
      % envsTrain and envsTest are cell arrays of environment lists
      
      % This function fills in all the necessary atom and bond contexts
      % into the models
      [atypes, atomContexts, bondContexts] = ...
         fillInContexts(modsTrain,envsTrain,modsTest,envsTest,includeAdhoc)
      % This function returns a Fitme object for the fits
      [ftrain ftest] = makeFitme(modsTrain,envsTrain,HLTrain, ...
         modsTest,envsTest,HLTest,includeAdhoc)
      % returns a rho(1,1) rho(1,2) rho(2,2) for the bond
      res = analyzeBond(mod,ienv,atom1,atom2)
      % Convert atom type to Z (used to determine ndim)
      res = atypeToZtype(atype)
      % 
   end % static methods
   methods
      function obj = Context(ndataIn,ztype1,ztype2)
         % preconfigures the storage assuming addModel will be called
         % n mod times
         if (nargin < 3) % atom context
            obj.isBondContext = 0;
            switch ztype1
               case 1
                  obj.ndim = 3;
               case 6
                  obj.ndim = 12;
            end
         else % bond context
            obj.isBondContext = 1;
            switch ztype1+ztype2
               case 2 % H2
                  obj.ndim = 3;
               case 7 % CH
                  obj.ndim = 12;
               case 12 % CC
                  obj.ndim = 7 * 3;
            end
         end
         obj.data = zeros(ndataIn,obj.ndim);
         obj.idata = 0;
      end
      function extractFeatures(obj)
         % performs feature extraction on data that has been added
         [dataz,obj.mu,obj.sigma] = zscore(obj.data);
         [obj.coeff,obj.score,obj.latent] = princomp(dataz);
      end
      function plotLatent(obj,figNum)
         figure(figNum);
         subplot(1,2,1);
         plot(obj.latent,'rx');
         title('latent');
         subplot(1,2,2);
         plot(cumsum(obj.latent)./sum(obj.latent),'bx');
         title('% variance explained');
      end
   end % methods
end

