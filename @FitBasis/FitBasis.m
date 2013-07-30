classdef FitBasis < handle
   % Adjust the definition of a low-level basis to agree with high-level
   properties(SetAccess = private)
      LLtplFile  % string with name of the template file
      HLtplFile  % string with name of the template file
      dataDir  % holds the tplFile and the cached data
      natom    % number of atoms (parsed from template file)
      HLnpar     % number of parameters in the template file
      LLnpar
      geomRandomSpecs % {} of same length as geomPars with structures
                      % that have the fields: name, low, high
      geomPars % {} of parameters in the tplFile that specify geometry
      extraPars % current values for the non-geometric parameters
      envs      % {} of environments
      HLbasis   % string specifying the HL basis set
      HLvalid   % false if HL data needs to be regenerated
      LLvalid   % false if LL data needs to be regenerated
      HL        % HL results
      LL        % LL results
      spin      % used in call to err
   end
   properties
      cacheHL   % useCache passed to Fragment for HL calc
      cacheLL   % as above for LL calc
   end      
   methods
      function obj = FitBasis(dataDir, HLtplFile, LLtplFile)
         obj.dataDir = dataDir;
         obj.HLtplFile = HLtplFile;
         obj.LLtplFile = LLtplFile;
         obj.geomRandomSpecs = {};
         obj.geomPars = {};
         obj.HLvalid = false;
         obj.LLvalid = false;
         obj.cacheHL = true;
         obj.cacheLL = false;
         LLtemplateText = fileread([dataDir,filesep,...
             LLtplFile,'.tpl']);
         HLtemplateText = fileread([dataDir,filesep,...
             HLtplFile,'.tpl']);
         obj.natom = size( strfind(LLtemplateText, 'ATOM'), 2);
         obj.HLnpar = size( strfind(HLtemplateText, 'PAR'), 2);
         obj.LLnpar = size( strfind(LLtemplateText, 'PAR'), 2);
         obj.spin = 1;
      end
      function setRandomGeom(obj,name,low,high, ntimes)
         % name = string for parameter, used only for display
         % low  = min for the random range
         % high = max for the random range
         % ntimes = add this ntimes (defaults to 1)
         if (nargin < 5)
            ntimes = 1;
         end
         t1.name = name;
         t1.low = low;
         t1.high = high;
         for i = 1:ntimes
            obj.geomRandomSpecs{end+1} = t1;
         end
         obj.HLvalid = false;
         obj.LLvalid = false;
      end
      function setExtraPars(obj,extraPars)
         obj.extraPars = extraPars;
         obj.LLvalid = false;
      end
      function res = getGeomPars(obj)
          res = obj.geomPars;
      end
      function addEnv(obj,env)
         obj.envs{end+1} = env;
         obj.HLvalid = false;
         obj.LLvalid = false;
      end
      function setHLBasis(obj, HLBasis)
         obj.HLbasis = HLBasis;
         obj.HLvalid = false;
      end
      function res = HLEnergies(obj, spin)
         if (~obj.HLvalid)
            obj.generateHLData(spin);
         end
         res = [];
         for i = 1:length(obj.HL)
            res = [res obj.extractEnergies(obj.HL{i})];
         end
      end
      function res = LLEnergies(obj, spin)
         if (~obj.LLvalid)
            obj.generateLLData(spin);
         end
         res = [];
         for i = 1:length(obj.LL)
            res = [res obj.extractEnergies(obj.LL{i})];
         end
      end
      function res = extractEnergies(obj,frag)
         res = [frag.EKE];
         for iatom = 1:frag.natom
            res = [res, frag.Een(iatom)];
         end
         res = [res, frag.E2];
      end
      function setSpin(obj,spin)
         obj.spin = spin;
      end
      function [res, hl, ll]  = err(obj,pars)
         obj.setExtraPars(pars);
         hl = obj.HLEnergies(obj.spin);
         ll = obj.LLEnergies(obj.spin);
         res = hl - ll;
         disp(['err(pars=',num2str(pars(:)'),') = ',num2str(norm(res))]);
      end
   end  
end

