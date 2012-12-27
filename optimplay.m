clear classes

mtrain = MSet;
%dfile = 'datasets/h2Dat.mat';
%mtrain.addData(dfile, 1:2:7, 1:10,1,796);
dfile = 'datasets/ch4rDat.mat';
mtrain.addData(dfile, 1:2, 1:4,1,791);
mtest = MSet;
mtest.addData(dfile, 5:6, 5:8,1,791);

fact  = MFactory;
%fact.setPolicies('h2fits');
fact.setPolicies('hybridslater1');
fact.makeMixInfo(mtrain.atomTypes);
f1 = fact.makeFitme(mtrain);
f1.silent = 1;
ftest = fact.makeFitme(mtest);
ftest.silent = 1;

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

monitor = OptMonitor(ftest);
monitor.plotNum = 300;
start = f1.getPars;
options = optimset('DiffMinChange',1.0e-4,'TolFun',1.0e-10, ...
   'TolX',3.0e-10,'MaxFunEvals',1e5,'Display','iter',...
   'OutputFcn',@monitor.toCall);
[pt,resnorm,residual,exitflag,output,lambda,jacobian] = ...
   lsqnonlin(@f1.err, start,lowLimits,highLimits,options);
monitor.storeoutput(pt,resnorm,residual,exitflag,output,lambda,jacobian);

%% what would a test set do
clear all;
load('optplay1.mat');
dfile = 'datasets/ch4rDat.mat';
mtest = MSet;
mtest.addData(dfile, 5:8, 5:8,1,791);
fact  = MFactory;
%fact.setPolicies('h2fits');
fact.setPolicies('hybridslater1');
fact.makeMixInfo(mtest.atomTypes);
ftest = fact.makeFitme(mtest);
ftest.silent = 1;

niter = length(monitor.values);
for i =1:niter
   display([num2str(i),' of ',num2str(niter)]);
   pars = monitor.param{i};
   x(i) = monitor.values{i}.iteration;
   etrain1(i) = monitor.values{i}.resnorm;
   etrain(i) = norm(f1.err(pars));
   etest(i) = norm(ftest.err(pars));
end

%% plot of residual
clear classes;
load('optplay1.mat');
[y,x] = monitor.trainError;
figure(100)
plot(x,y,'bo');
[y,x] = monitor.testError;
hold on;
plot(x,y,'ro');
figure(200)
plot(-diff(y),'ro');

%% print things out
% clear all;
% load('optplay1.mat');
clc;
niter = length(monitor.values);
for i=2:niter
   vp = monitor.values{i-1};
   v = monitor.values{i};
  disp(['ITERATION ',num2str(i),' type ',monitor.state{i}]);
  r1 = sqrt(v.resnorm);
  r2 = norm(v.residual);
  disp(['norm(v.residual) ',num2str(r2),' sqrt(v.resnorm) ',num2str(r2),...
     ' norm(f1.err) ',num2str(etrain(i))]);
  disp(['vdiff ',num2str(v.resnorm-vp.resnorm),...
     'xdiff ',num2str(max(abs(monitor.param{i}-monitor.param{i-1})))]);
  v
end

