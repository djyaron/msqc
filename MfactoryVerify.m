%% Create original style 
clear classes
%load('datasets/ch4rDat.mat');
load('datasets/ethanerDat.mat');
mtrain = cell(0,0);
HLtrain = cell(0,0);
envsTrain = cell(0,0);
for i = 1:2
   mtrain{end+1} = Model3(LL{i,1},LL{i,1},LL{i,1});
   mtrain{end}.solveHF;
   HLtrain{end+1} = HL{i,1};
   envsTrain{1,end+1} = 1:4;
end
mtest = cell(0,0);
HLtest = cell(0,0);
envsTest = cell(0,0);
for i = 3:4
   mtest{end+1} = Model3(LL{i,1},LL{i,1},LL{i,1});
   mtest{end}.solveHF;
   HLtest{end+1} = HL{i,1};
   envsTest{1,end+1} = 1:4;
end
includeAdhoc = 1;
separateSP = 0;
include1s = 0;
hybrid = 1;
[forig ftest] = Context.makeFitme(mtrain,envsTrain,HLtrain, ...
   mtest,envsTest,HLtest,includeAdhoc,separateSP,include1s,hybrid);
save('debug1.mat','forig');

%%
clear classes;
load('debug1.mat');
ms = MSet;
%ms.addData('datasets/ch4rDat.mat',1:2,1:4,1,791);
ms.addData('datasets/ethanerDat.mat',1:2,1:4,1,791);
disp(['atom types ',num2str(ms.atomTypes)]);

m1 = MFactory;
m1.addPolicy('o','KE', 'f','scale', 'sp','sonly', 'i',1, 'c','q r bo');
m1.addPolicy('o','EN', 'f','scale', 'sp','sonly', 'i',1, 'c','q r bo');
m1.addPolicy('o','E2', 'f','scale', 'sp','sonly', 'i',1, 'c','q r bo');

m1.addPolicy('o','KE', 'f','scale', 'sp','hybrid', 'i',6, 'j',1, ...
   'c','r bo q');
m1.addPolicy('o','EN', 'f','scale', 'sp','hybrid', 'i',6, 'j',1, ...
   'c','r bo q');
m1.addPolicy('o','E2', 'f','scale', 'sp','hybrid', 'i',6, 'j',1, ...
   'c','r bo q');

m1.addPolicy('o','KE', 'f','scale', 'sp','combine', 'i',6, 'c','q r bo');
m1.addPolicy('o','KE', 'f','const', 'i',6, 'sp','combine');
m1.addPolicy('o','EN', 'f','scale', 'sp','combine', 'i',6, 'c','q r bo');
m1.addPolicy('o','EN', 'f','const', 'i',6, 'sp','combine');
m1.addPolicy('o','E2', 'f','scale', 'sp','combine', 'i',6, 'c','q r bo');

m1.addPolicy('o','KE', 'f','scale', 'sp','hybrid', 'i',6, 'j',6, ...
   'c','r bo q');
m1.addPolicy('o','EN', 'f','scale', 'sp','hybrid', 'i',6, 'j',6, ...
   'c','r bo q');
m1.addPolicy('o','E2', 'f','scale', 'sp','hybrid', 'i',6, 'j',6, ...
   'c','r bo q');

m1.addPolicy('o','E2', 'f','scale', 'sp','sonly', 'i',1, 'j',1, ...
   'c','r');

m1.makeMixInfo(ms.atomTypes);
m1.printMixInfo;

[f1,c1] = m1.makeFitme(ms);

%%
for i = 1:length(f1.mixers)
   if (length(f1.mixers{i}.fixed) > 1)
      f1.mixers{i}.fixed = zeros(size(f1.mixers{i}.fixed));
   end
end

for i = 1:length(forig.mixers)
   if (length(forig.mixers{i}.fixed) > 2)
      for j= 1:3
         forig.mixers{i}.fixed(2:4) = 0;
      end
   else
      forig.mixers{i}.fixed(2) = 0;
   end
end


%%
e1 = f1.err(f1.getPars);
f1.printEDetails;

forig.includeEtot = 1;
forig.silent = 1;
forig.plot = 0;
eorig = forig.err(forig.getPars);
forig.printEDetails;
disp(['starting diff ',num2str(max(abs(e1-eorig)))]);
%%
dis = zeros(f1.npar,1);
for ipar = 1:f1.npar
   xsave = f1.getPars;
   x = xsave;
   x(ipar) = x(ipar) + 0.1;
   e1 = f1.err(x);
   f1.setPars(xsave);
   xsave = forig.getPars;
   x = xsave;
   x(ipar) = x(ipar) + 0.1;
   eorig = forig.err(x);
   forig.setPars(xsave);
   dis(ipar) = max(abs(eorig-e1));
   if (dis(ipar) > 0.001)
      xsave = f1.getPars;
      x = xsave;
      x(ipar) = x(ipar) + 0.1;
      f1.setPars(x);
      f1.models{1}.printMixers;
      f1.setPars(xsave);
      xsave = forig.getPars;
      x = xsave;
      x(ipar) = x(ipar) + 0.1;
      forig.setPars(x);
      forig.models{1}.printMixers;
      forig.setPars(xsave);
   end
   disp(['ipar = ',num2str(ipar),' diff ',num2str(dis(ipar))]);
end
%%
f1.silent = 0;
lowLimits = zeros(f1.npar,1);
highLimits = lowLimits;
i1 = 1;
for imix = 1:length(f1.mixers)
   mix = f1.mixers{imix};
   if (strcmp(mix.funcType,'scale'))
      lowLimits(i1) = -1.0;
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
