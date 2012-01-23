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
config.par = [par 1.05 1.05 1.05 1.05 1.05];
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
%% Fit to reg
mod = Model2(reg,nar,dif);
f1 = Fitme;
f1.addFrag(mod,reg);
pt{1} = zeros(1,16);
for i=1:4
   disp(['updating density for p= ',num2str(pt{i})]);
   f1.updateDensity(pt{i});
   err{i,1} = f1.err(pt{i});
   disp(['rms err = ',num2str(sqrt(err{i,1}*err{i,1}'))]);
   pt{i+1} = lsqnonlin(@f1.err, pt{i});
   err{i+1,2} = f1.err(pt{i+1});
   disp(['rms err approx = ',num2str(sqrt(err{i+1,2}*err{i+1,2}'))]);
end




