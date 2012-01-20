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

%mod = Model2(reg,nar,dif,reg);
mod = Model2(reg,nar,dif,reg);
mod.par = zeros(1,16);

f1 = Fitme;
f1.addFrag(mod,reg);
p1 = zeros(1,16);
f1.updateDensity(p1);
%%
 p2 = lsqnonlin(@f1.err, p1);
 %%
 f1.updateDensity(p2);
 p3 = lsqnonlin(@f1.err, p2);




