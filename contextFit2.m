function [err,pt,testErr] = contextFit2(f1,ftest,mixNumber,parNumber,...
   remove,maxIter,psc)

if (nargin < 5)
   remove = 0;
end
if (nargin < 6)
   maxIter = 500;
end
if (nargin < 7)
   psc = 0;
end
if (isempty(f1))
   load('f1temp.mat','f1','ftest');
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
ftest.plot = 0;
% dataDir = [topDir,filePre,'/cfit/'];
% if (exist(dataDir,'dir') ~= 7)
%    status = mkdir(dataDir);
% end

%   diary([dataDir,'out.diary']);
%   diary on;
tic

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
if (psc == 0)
   options = optimset('DiffMinChange',1.0e-4,'TolFun',1.0e-3, ...
      'TolX',3.0e-3,'MaxFunEvals',maxIter);
   [pt,resnorm,residual,exitflag,output,lambda,jacobian] = ...
      lsqnonlin(@f1.err, start,lowLimits,highLimits,options);
   err = sqrt(residual*residual'/length(residual))*627.509;
elseif (psc == 1)
   options = {'Display',0,'FunTol',1.0e-3,'XTol',1.0e-3,'Lambda',1.0e-6};
   [pt,ssq, CNT, Res, XY] = LMFnlsq(@f1.err,start',options{:});
   pt = pt';
   err = sqrt(ssq/f1.ndata)*627.509;
else
   ista(@f1.err,start);
end
testres = ftest.err(pt);
testErr = sqrt(testres*testres'/length(testres))*627.509;
%pt
%resnorm
%f1.printMixers;

end