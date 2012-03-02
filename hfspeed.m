% Scripts for working on speeding up the Hartree Fock algorithm
%   In Model3, copied 

%% load some data
clear classes
root='C:\dave\apoly\msqc\';
load('ethane4/env2.mat');
config = Fragment.defaultConfig();
config.method = 'MP2';
config.par = [1.54 1.12 60];
config.template = 'ethane1';
config.basisSet = 'STO-3G';
frag = Fragment([root,'ethane4mp2'], config);
for ienv = 1:2
   frag.addEnv(env{ienv});
end


m1 = Model3(frag,frag,frag);
m2 = Model3(frag,frag,frag);

tic;
[a.orb,a.Eorb,a.Ehf] = m1.hartreeFockSlow(1,1e-8,100,50);
t1 = toc;
tic;
[b.orb,b.Eorb,b.Ehf] = m2.hartreeFock(1,1e-8,100,50);
t2 = toc;

disp(['time 1 = ',num2str(t1),' time 2 = ',num2str(t2), ...
   ' diff Eorb ',num2str(max(max(abs(a.Eorb-b.Eorb)))), ...
   ' diff Ehf ',num2str(abs(a.Ehf-b.Ehf))]);

