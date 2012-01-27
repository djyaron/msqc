%% Get some test fragments
clear all;
root = 'c:\dave\apoly\msqc\';
load('ethane4/env2.mat');
config = Fragment.defaultConfig();
par = [1.54 1.12 60];
config.par = par;
config.template = 'ethane1';
config.basisSet = 'STO-3G';
reg = Fragment([root,'ethane4'], config);
config.template = 'ethane1-gen';
config.basisSet = 'GEN';
config.par = [par 0.9 0.9 0.9 0.9 0.9];
nar = Fragment([root,'ethane4'], config);
config.par = [par 1.1 1.1 1.1 1.1 1.1];
dif = Fragment([root,'ethane4'], config);
for ienv = 1:20 %size(env,2)
   reg.addEnv(env{ienv});
   nar.addEnv(env{ienv});
   dif.addEnv(env{ienv});
end
save('ethane4\m2verify.mat');
%%
clear classes;
load('ethane4\m2verify.mat');
mod = Model2(reg,nar,dif);
mod.par = rand(1,16);
%% test of the new H1
h1diff = max(max(abs(mod.H1check(0) - mod.H1(0))))
%% test of range on partitionE1 (to avoid the sum(sum()) thing).
mod.updateDensity(mod.par);
t1 = sum(sum( mod.partitionE1(0) ) );
arange = cell(1,1);
arange{1,1}= 1:mod.nbasis;
t2 = mod.partitionE1(0,mod.H1(0),arange);
tdiff = t1-t2
%% Test of fitting: Should get zero error for this set-up
mod = Model2(reg,reg,reg);
f1 = Fitme;
f1.addFrag(mod,reg);
p1 = rand(1,16);
f1.updateDensity(p1);
t3 = max(abs(f1.err(p1)));
%% Fit to reg (failed because it does a fit to 16 params from start)
clear classes;
load('ethane4\m2verify.mat');
mod = Model2(reg,nar,dif);
f1 = Fitme;
f1.addFrag(mod,reg);
pt{1} = zeros(1,16);
for i=1:4
   disp(['updating density for p= ',num2str(pt{i})]);
   f1.updateDensity(pt{i});
   err{i,2} = f1.err(pt{i});
   disp(['rms err = ',num2str(sqrt(err{i,2}*err{i,2}'))]);
   pt{i+1} = lsqnonlin(@f1.err, pt{i});
   err{i+1,1} = f1.err(pt{i+1});
   disp(['rms err approx = ',num2str(sqrt(err{i+1,1}*err{i+1,1}'))]);
end
%% Fit by ramping up number of parameters
clear classes;
load('ethane4\m2verify.mat');
mod = Model2(reg,nar,dif);
f1 = Fitme;
f1.addFrag(mod,reg);
f1.exactDensity = 1;
% start with 4 parameter fit
mod.sepKE = 0;
mod.sepSP = 0;
pt{1} = zeros(1,mod.npar);
scaler = ones(1,mod.npar);
for i=1:1
   disp(['updating density for p= ',num2str(pt{i})]);
   f1.updateDensity(pt{i});
   err{i,2} = f1.err(pt{i});
   disp(['rms err = ',num2str(sqrt(err{i,2}*err{i,2}'))]);
%   options = optimset('MaxIter',1);
%   pt{i+1} = lsqnonlin(@f1.err, pt{i},[],[],options);
   %bounds = 0.2 * i *scaler;
   pt{i+1} = lsqnonlin(@f1.err, pt{i});% ,-bounds,bounds);  
   err{i+1,1} = f1.err(pt{i+1});
   disp(['rms err approx = ',num2str(sqrt(err{i+1,1}*err{i+1,1}'))]);
end
%% Plot error versus p
clear classes;
load('ethane4\m2verify.mat');
mod = Model2(reg,nar,dif);
ic = 0;
for ps1 = -1:0.2:1
   ic = ic+1;
   pars = ps1 * ones(1,16);
   mod.setPar(pars);
   mod.solveHF();
   p5(ic) = ps1;
   for ienv = 0:mod.nenv
      ke5(ienv+1,ic) = sum(sum(mod.partitionE1(ienv,mod.KE)));
   end
end
%%
for ienv = 0:reg.nenv
   ke0(ienv+1,1) = sum(sum(reg.partitionE1(ienv,reg.KE)));
end
%%
for i=1:reg.nenv+1
   hold on;
   plot(p5,ke5(i,:)+i,'b-');
   hold on;
   plot(0,ke0(i,1)+i,'rx');
end

%%
clear classes;
load('ethane4\m2verify.mat');
mod = Model2(reg,nar,dif);
f1 = Fitme;
f1.addFrag(mod,reg);
f1.exactDensity = 1;
% start with 4 parameter fit
mod.sepKE = 0;
mod.sepSP = 0;
pt{1} = zeros(1,mod.npar);
f1.updateDensity(pt{1});
err{1} = f1.err(pt{1});
disp([' starting rms err = ',num2str(sqrt(err{1}*err{1}'))]);
pt{2} = lsqnonlin(@f1.err, pt{1});
err{2} = f1.err(pt{2});
disp(['KE=0 SP=0 rms err = ',num2str(sqrt(err{2}*err{2}'))]);
%%
mod.sepKE = 1;
pstart = mod.mapPar( mod.par );
pt{3} = lsqnonlin(@f1.err, pstart);
err{3} = f1.err(pt{3});
disp(['KE=1 SP=0 rms err = ',num2str(sqrt(err{3}*err{3}'))]);
%% reset back to 4 parameter fit, and then try sepSP = 1
mod.sepKE = 0;
mod.setPar(pt{2});
mod.sepSP = 1;
pstart = mod.mapPar( mod.par );
pt{4} = lsqnonlin(@f1.err, pstart);
err{4} = f1.err(pt{4});
disp(['KE=0 SP=1 rms err = ',num2str(sqrt(err{4}*err{4}'))]);
%% full 16 parameter fit
mod.setPar(pt{4}); % set parameters based on above fit
mod.sepSP = 1; % adjust new assumptions
mod.sepKE = 1;
pstart = mod.mapPar( mod.par );
pt{5} = lsqnonlin(@f1.err, pstart);
err{5} = f1.err(pt{5});
disp(['KE=1 SP=1 rsms err = ',num2str(sqrt(err{4}*err{4}'))]);


