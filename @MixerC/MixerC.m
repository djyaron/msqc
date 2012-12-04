classdef MixerC < handle
   properties
      par     % (1,npar) current parameters
      desc    % string description
      fixed   % (1,npar) 0 if parameter should be fit, 1 if fixed
      hybrid  % 0 = no hybridization, 1 = sigma mods, 2 = pi mods
      isDiag  % 1 if diagonal
      funcType % scale const interp
      bonded  % true = bonded only  false = nonbonded only
      context % string specifying contexts
      cset    % pointer to object with method context(imodel,imixer)
      index   % used and managed externally
   end

   methods
      function obj = MixerC(par,funcType,hybrid)
         if (nargin < 1)
            par = [0];
         end
         if (nargin < 2)
            funcType = 'scale';
         end
         if (nargin < 3)
            hybrid = 0;
         end
         obj.par = par;
         obj.desc = '';
         obj.fixed = ones(size(par));
         obj.hybrid = hybrid;
         obj.bonded = true;
         obj.funcType = validatestring(funcType,{'scale','const','interp'});
         obj.context = [];
      end
      function res = deepCopy(obj)
         res = Mixer(obj.par,obj.funcType,obj.hybrid);
         res.fixed = obj.fixed;
         res.desc = obj.desc;
      end
      function res = npar(obj)
         res = sum(obj.fixed == 0);
      end
      function res = setPars(obj, pars)
         ic = find(obj.fixed == 0);
         if (length(ic) > 0)
            obj.par(ic) = pars;
         end
      end
      function res = getPars(obj)
         ic = find(obj.fixed == 0);
         res = obj.par(ic);
      end
      function res = mixFunction(obj,x,v0,v1,v2, model, ii, jj)
         if (obj.hybrid)
            % mixFunctionHybrid rotates to hybrid orbitals
            % and then call mixFunctionNormal on these rotated orbitals
            res = obj.mixFunctionHybrid(x,v0,v1,v2, model, ii, jj);
         else
            res = obj.mixFunctionNormal(x,v0,v1,v2);
         end
      end
      function res = mixFunctionNormal(obj,x,v0,v1,v2)
         switch obj.funcType
            case 'scale'
               res = (1.0 + x) * v0;
            case 'const'
               if (length(size(v0))==4)
                  res = v0;
                  for i = 1:size(v0,1)
                     res(i,i,i,i) = res(i,i,i,i) + x;
                  end
               else
                  res = v0 + x * eye(size(v0));
               end
            case 'interp'
               res = ((1.0-x)/2.0) * v1 + ((1.0+x)/2.0) * v2;
            otherwise
               error(['Mixer: mixFunctionNormal: unknown functype ',...
                  num2str(obj.funcType)]);
         end
      end
      function res = mix(obj, v0, v1, v2, model, ii, jj, ienv)
         iatom = model.basisAtom(ii(1));
         jatom = model.basisAtom(jj(1));
         bond = model.isBonded(iatom,jatom);
         if ((iatom ~= jatom) && (obj.bonded ~= bond))
            res = v0;
            return
         end
         if (isempty(obj.context))
            x = obj.par(1);
         else
            contextVars = ...
               obj.cset.getContext(model,obj,iatom,jatom,ienv);
            contextPars = obj.par(2:end);
            x = obj.par(1) + sum(contextPars(:) .* contextVars(:));
         end
         res = obj.mixFunction(x,v0,v1,v2, model, ii, jj);
      end
      function res = toString(obj)
         htypes = {'','hyb-sig','hyb-pi'};
         res = [obj.desc,' ',obj.funcType,' ',htypes{obj.hybrid+1}];
         if (~obj.bonded)
            res = [res,' nonbond'];
         end
         res = [res,' ', obj.context];
         for i = 1:size(obj.par,2)
            res = [res,' ',num2str(obj.par(i))];
            if (obj.fixed(i))
               res = [res,'*'];
            end
         end
      end
      function print(obj)
         disp(obj.toString);
      end
   end
end

