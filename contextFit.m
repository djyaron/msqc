function [err,pt] = contextFit(f1,mixNumber,parNumber,remove, maxIter)

if (nargin < 4)
   remove = 0;
end
if (nargin < 5)
   maxIter = 500;
end
if (isempty(f1))
   load('f1temp.mat','f1');
end
if (mixNumber > 0)
   if (remove == 0)
      f1.mixers{mixNumber}.fixed(parNumber) = 0;
   else
      f1.mixers{mixNumber}.fixed(parNumber) = 1;
      f1.mixers{mixNumber}.par(parNumber) = 0.0;
   end
end

f1.plot = 0; % showPlots;
f1.parallel = 0;
% dataDir = [topDir,filePre,'/cfit/'];
% if (exist(dataDir,'dir') ~= 7)
%    status = mkdir(dataDir);
% end

%   diary([dataDir,'out.diary']);
%   diary on;
tic
options = optimset('DiffMinChange',1.0e-4,'TolFun',1.0e-3, ...
   'TolX',3.0e-3,'MaxFunEvals',maxIter);
lowLimits = zeros(f1.npar,1);
highLimits = lowLimits;
i1 = 1;
for imix = 1:length(f1.mixers)
   mix = f1.mixers{imix};
   if ((mix.funcType == 2)||(mix.funcType == 3))
      lowLimits(i1) = 0.0;
      highLimits(i1) = 10.0;
      i1 = i1+1;
      for i2 = 2:mix.npar
         lowLimits(i1) = -inf;
         highLimits(i1) = inf;
         i1 = i1+1;
      end
   else
      for i2 = 1:mix.npar
         lowLimits(i1) = -inf;
         highLimits(i1) = inf;
         i1 = i1+1;
      end
   end
end
start = f1.getPars;
[pt,resnorm,residual,exitflag,output,lambda,jacobian] = ...
   lsqnonlin(@f1.err, start,lowLimits,highLimits,options);
err = sqrt(residual*residual'/length(residual))*627.509;
%pt
%resnorm
%f1.printMixers;

end