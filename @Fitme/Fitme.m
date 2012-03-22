classdef Fitme < handle
   properties
      models  % {1,nmodels}  cell array of models
      HLs     % {1,nmodels}  cell array to fit to
      HLKE    % {1,nmodels}(1,nenv) KE energy
      HLEN    % {1,nmodels}(natom,nenv) electron-nuclear interaction
      mixers  % {1,nmixer}   cell
      
      envs      % {1,nmodels} list of environments to include in fit
      includeKE % include kinetic energy in fit
      includeEN % {1,Z} include elec-nuc operators for element Z
      
      parHF   % Last parameters for which HF was solved
      epsDensity % re-evaluate density matrix if par change > eps
   end
   methods
      function res = Fitme
         res.models = cell(0,0);
         res.HLs    = cell(0,0);
         res.HLKE   = cell(0,0);
         res.epsDensity = 0.0;
         res.includeKE = 1;
         res.includeEN = zeros(1,6);
         res.parHF = [];
      end
      function addMixer(obj, mix)
         add = 1;
         for i=1:size(obj.mixers,2)
            if (mix == obj.mixers{1,i})
               add = 0;
               break;
            end
         end
         if (add == 1)
            obj.mixers{1,end+1} = mix;
         end
      end
      function addMixers(obj,mixersIn)
         % mixes is a cell array (1,:) of mixers
         for i=1:size(mixersIn,2)
            obj.addMixer(mixersIn{i});
         end
      end
      function addFrag(obj,model,HL)
         obj.models{1,end+1} = model;
         obj.addMixers(model.mixers);
         obj.HLs{1,end+1} = HL;
      end
      function setEnvs(obj,envsIn)
         % currently assumes same environments for every frag/model
         obj.envs = cell(1,obj.nmodels);
         for i=1:obj.nmodels
            obj.envs{1,i} = envsIn;
         end
         obj.HLKE = cell(0,0);
         obj.HLEN = cell(0,0);
         for imod = 1:obj.nmodels
            envs1 = obj.envs{1,i};
            HL = obj.HLs{imod};
            obj.HLKE{1,end+1} = HL.EKE(envs1);
            nenv = size(envs1,2);
            en = zeros(HL.natom, nenv);
            for iatom = 1:HL.natom
               en(iatom,:) = HL.Een(iatom,envs1);
            end
            obj.HLEN{1,end+1} = en;
         end
         obj.parHF = [];
      end
      function res = nmodels(obj)
         res = size(obj.models,2);
      end
      function res = npar(obj)
         res = 0;
         for i=1:size(obj.mixers,2)
            res = res + obj.mixers{i}.npar;
         end
      end
      function res = getPars(obj)
         res = zeros(1,obj.npar);
         ic = 1;
         for i = 1:size(obj.mixers,2)
            mtemp = obj.mixers{1,i};
            np = mtemp.npar;
            if (np > 0)
               res(ic:(ic+np-1)) = mtemp.getPars;
            end
            ic = ic + np;
         end
      end
      function setPars(obj,par)
         % sets parameters, and updates densities
         if (size(par,2) ~= obj.npar)
            error(['Fitme.SetPars called with ',num2str(size(par,2)), ...
               ' parameters when ',num2str(obj.npar),' are needed ']);
         end
         ic = 1;
         for i = 1:size(obj.mixers,2)
            mtemp = obj.mixers{1,i};
            np = mtemp.npar;
            if (np > 0)
               mtemp.setPars( par(ic:(ic+np-1)));
            end
            ic = ic + np;
         end
      end
      function updateDensity(obj)
         par = obj.getPars;
         if (size(obj.parHF,1) == 0)
            dpar = 1e10;
         else
            dpar = max(abs(obj.parHF-par));
         end
         if (dpar > obj.epsDensity)
            disp(['solving for density matrices']);
            for imod = 1:obj.nmodels
               obj.models{imod}.solveHF(obj.envs{1,imod});
            end
            obj.parHF = par;
         end
      end
      function res = ndata(obj)
         ic = 0;
         for imod = 1:obj.nmodels
            if (obj.includeKE == 1)
               ic = ic + size(obj.HLKE{1,imod},2);
            end
            for iatom = 1:obj.HLs{imod}.natom
               if (obj.includeEN( obj.HLs{imod}.Z(iatom) ))
                  ic = ic + size(obj.HLEN{imod}(iatom,:),2);
               end
            end
         end
         res = ic;
      end
      function res = err(obj,par)
         flip = 0; % to handle fit routines that pass row or column
         if (size(par,1)>size(par,2))
            par = par';
            flip = 1;
         end
         disp(['Fitme.err called with par = ',num2str(par)]);
         obj.setPars(par);
         obj.updateDensity();

         ic = 1;
         ndat = obj.ndata;
         res = zeros(1,ndat);
         for imod = 1:obj.nmodels
            if (obj.includeKE == 1)
               t1 = obj.HLKE{1,imod} - ...
                  obj.models{imod}.EKE(obj.envs{1,imod});
               n = size(t1,2);
               res(1,ic:(ic+n-1))= t1;
               ic = ic + n;
            end
            for iatom = 1:obj.HLs{imod}.natom
               if (obj.includeEN( obj.HLs{imod}.Z(iatom) ))
                  t1 = obj.HLEN{imod}(iatom,:) - ...
                     obj.models{imod}.Een(iatom,obj.envs{1,imod});
                  n = size(t1,2);
                  res(1,ic:(ic+n-1)) = t1;
                  ic = ic + n;
               end
            end
         end
         disp(['RMS err/ndata = ',num2str(sqrt(res*res')/ndat)]);
         if (flip == 1)
            res = res';
         end
      end
   end   
end

