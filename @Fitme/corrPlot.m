function [LL,HL] = corrPlot(obj, par, oper, figNum, color)
% oper = 0 (KE)
% oper = # (electron nuclear of nth nucleus

if (nargin < 3)
   oper = 0;
end
if (nargin < 4)
   figNum = 800 + oper;
end
if (nargin < 5)
   color = 'b';
end

% update models
for imod = 1:obj.nmodels
   obj.models{imod}.setPar(par);
end
if (size(obj.parHF,1) == 0)
   dpar = 1;
else
   dpar = max(abs(obj.parHF-par));
end
if (obj.exactDensity && (dpar > 1.0e-4))
   disp(['solving for density matrices']);
   for imod = 1:obj.nmodels
      obj.models{imod}.solveHF();
   end
   obj.parHF = par;
end

% Do sum over all orbitals
sumRange = cell(1,1);
ic = 0;
for imod = 1:obj.nmodels
   sumRange{1,1} = 1:obj.models{imod}.nbasis;
   if (oper == 0)
      for ienv = 0:obj.models{imod}.nenv
         ic = ic + 1;
         LL(ic) = obj.models{imod}.partitionE1(ienv , ...
            obj.models{imod}.KE, sumRange);
         HL(ic) = obj.HLKE{1,imod}(ienv+1);
      end
   end
   if (oper > 0)
      iatom = oper;
      for ienv = 0:obj.models{imod}.nenv
         ic = ic+1;
         LL(ic) = obj.models{imod}.partitionE1(ienv, ...
            obj.models{imod}.H1en(iatom), sumRange);
         HL(ic) = obj.HLEN{1,imod}(iatom,ienv+1);
      end
   end
end
figure(figNum);
plot(HL,LL,[color,'.']);

if (oper == 0)
   title('Kinetic energy');
else
   title(['Electron Nuc to atom ',num2str(iatom)]);
end

