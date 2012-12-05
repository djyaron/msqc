classdef MFactory < handle
    % Builds models in a multistage process
    % 1) Make an empty MFactory object and call addPolicy to add
    %    rules used to construct mixers. The rules specify operator
    %    sp (hybrid, sonly, ponly) and atom types (iatom for diagonal,
    %    iatom and jatom for off diagonal). The input to this function is
    %    varargin, with 'key','value' pairs. the possible keys and their
    %    values are specified in the verify* methods called at the start of
    %    the addPolicy method.
    % 2) Call MFactory.makeMixInfo(atomtypes) to create mixInfo{} which
    %    holds all mixers that may be needed for the given set of
    %    atomtypes. mixInfo includes both the mixer and the information
    %    regarding how that mixer should be added to models. (Note that
    %    MSet has a method "atomTypes" that returns all atom times in a
    %    given set of models. 
    % 3) Call MFactory.makeFitme(MSet) to creat a fitme object that
    %    contains all of the models in MSet. makeFitme: adds the models,
    %    call Fitme.setEnv() to calculate all of the highlevel results
    %    against which to fit, and creates a CSet object that has all of
    %    the needed contexts. (The context variables depend on both mixer
    %    and model, which is why this is handled here.)
    % The motivation for this three stage building is that stage 1 handles
    % the decisions (how is each operator being handled for each atom
    % type). In some cases, a decision may need to generate many mixers, so
    % it would not be good to handle these decisions through specific
    % mixers.
    % Stage 2 implements those decisions for a particular set of atom
    % types. This includes making the mixers, and since the mixers hold the
    % paramters, the mixers are actually functions of both the policies and
    % the training data. At the end of stage 2, you essentially have an
    % object that you can train on any set of data.
    % Stage 3 applies the current state of mixers to a set of models to
    % make a fitme object. To begin, you would make a fitme object based on
    % the training data, and then do a fit of the parameters. Once these
    % parameters are established, additional calls to makefitme can be used
    % to apply those trained parameters to any set of models.
   properties
      policy      % cell array of policies
      atomTypes   % atom types used to construct mixInfo
      mixInfo     % cell array of structures specifying mixers
      mixer       % array of unique mixers
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
            {'jatom','j'},'sp',{'context','contexts','c'},{'nonbond','nb'}};
         t1 = validateInput(varargin,validParameters);
         t2.oper = validatestring(t1.oper,{'KE','EN','E2','*'});
         t2.func = validatestring(t1.func,{'const','scale','interp'});
         t2.sp   = validatestring(t1.sp, ...
          {'core','sonly','ponly','hybrid','combine','separate','slater'});
         t2.iatom = t1.iatom;
         if (isfield(t1,'jatom'))
            t2.jatom = t1.jatom;
            if (isfield(t1,'nonbond'))
               t2.nonbond = true;
            else
               t2.nonbond = false;
            end
            % Need to flip these due to the way we avoid adding mixers
            % twice in makeMixInfo
            if (t2.iatom < t2.jatom)
                temp1 = t2.iatom;
                t2.iatom = t2.jatom;
                t2.iatom = temp1;
            end
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
