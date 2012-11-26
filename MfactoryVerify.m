%%
clear classes

ms = MSet;
ms.addData('datasets/ch4rDat.mat',1:10,1:2:20,1,791);
disp(['atom types ',num2str(ms.atomTypes)]);

m1 = MFactory;
% Constants
m1.addPolicy('o','KE', 'f','const', 'i',1, 'sp','sonly');
m1.addPolicy('o','KE', 'f','const', 'i',6, 'sp','sonly');
m1.addPolicy('o','KE', 'f','const', 'i',6, 'sp','ponly');

m1.addPolicy('o','EN', 'f','const', 'i',1, 'sp','sonly');
m1.addPolicy('o','EN', 'f','const', 'i',6, 'sp','sonly');
m1.addPolicy('o','EN', 'f','const', 'i',6, 'sp','ponly');

% diagonal
m1.addPolicy('o','KE', 'f','scale', 'sp','combine', 'i','*', 'c','r q bo');
m1.addPolicy('o','EN', 'f','scale', 'sp','combine', 'i','*', 'c','r q bo');
m1.addPolicy('o','E2', 'f','scale', 'sp','combine', 'i','*', 'c','r q bo');

% bonded
m1.addPolicy('o','KE', 'f','scale', 'sp','hybrid', 'i','*', 'j','*', ...
   'c','r q bo');
m1.addPolicy('o','EN', 'f','scale', 'sp','hybrid', 'i','*', 'j','*', ...
   'c','r q bo');
m1.addPolicy('o','E2', 'f','scale', 'sp','hybrid', 'i','*', 'j','*', ...
   'c','r q bo');

m1.makeMixInfo(ms.atomTypes);
m1.printMixInfo;

[f1,c1] = m1.makeFitme(ms);

e1 = f1.err(f1.getPars);
f1.printEDetails;

lowLimits = zeros(f1.npar,1);
highLimits = lowLimits;
i1 = 1;
for imix = 1:length(f1.mixers)
   mix = f1.mixers{imix};
   if (strcmp(mix.funcType,'scale'))
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
maxIter = 500;
options = optimset('DiffMinChange',1.0e-4,'TolFun',1.0e-3, ...
   'TolX',3.0e-3,'MaxFunEvals',maxIter);
[pt,resnorm,residual,exitflag,output,lambda,jacobian] = ...
   lsqnonlin(@f1.err, start,lowLimits,highLimits,options);
err = sqrt(residual*residual'/length(residual))*627.509;
