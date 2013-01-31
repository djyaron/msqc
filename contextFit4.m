function [err,pt, testErr, monitor] = contextFit4(f1,ftest,...
   maxIter,epsTest,updateContext)
if (nargin < 5)
   updateContext = 0;
end
lowLimits = zeros(f1.npar,1);
highLimits = lowLimits;
i1 = 1;
for imix = 1:length(f1.mixers)
   mix = f1.mixers{imix};
   if (strcmp(mix.funcType,'scale'))
      lowLimits(i1) = -1.0;
      highLimits(i1) = 2.0;
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

monitor = OptMonitor(f1,ftest,maxIter,epsTest);
monitor.plotNum = 99;
monitor.updateContext = updateContext;
start = f1.getPars;
options = optimset('DiffMinChange',1.0e-4,'TolFun',1.0e-10, ...
   'TolX',3.0e-10,'MaxFunEvals',1e5,'Display','iter',...
   'OutputFcn',@monitor.toCall);
[pt,resnorm,residual,exitflag,output,lambda,jacobian] = ...
   lsqnonlin(@f1.err, start,lowLimits,highLimits,options);
monitor.storeoutput(pt,resnorm,residual,exitflag,output,lambda,jacobian);

err = sqrt(residual*residual'/length(residual))*627.509;
if (~isempty(ftest))
   testres = ftest.err(ftest.getPars);
   testErr = sqrt(testres*testres'/length(testres))*627.509;
end
end