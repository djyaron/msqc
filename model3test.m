%%
clear classes;
load('datasets/ethaneDat.mat');
frag = LL{1,1};
fnar = LL{1,2};
fdif = LL{1,3};
fhl = HL{1,1};

save('temp.mat','frag','fnar','fdif','fhl');

%%
clear classes;
% Need to have frag, fnar and fdif stored in temp.mat
load('temp.mat');
m = Model3(frag,fnar,fdif);
m.addKEmodDiag(6,[1,2]);
m.addKEmodDiag(1,1);
m.addKEmodBonded(1,6,[1 2],[1 2]);
m.addKEmodBonded(6,6,[1 2],[1 2]);
mods = m.KEmods;
for i=1:size(mods,2)
   disp(['KE ',num2str(i),' ilist [',num2str(mods{i}.ilist), ...
      '] jlist [',num2str(mods{i}.jlist),'] ',mods{i}.mixer.desc]);
end
m.addENmodDiag(6,[1,2]);
m.addENmodDiag(1,1);
m.addENmodBonded(6,1,[1 2],[1 2]);
m.addENmodBonded(6,6,[1 2],[1 2]);
for iatom = 1:m.natom
   disp(['H1 EN for atom ',num2str(iatom)]);
   mods = m.ENmods{1,iatom};
   for i=1:size(mods,2)
      disp(['EN ',num2str(i),' ilist [',num2str(mods{i}.ilist), ...
         '] jlist [',num2str(mods{i}.jlist),'] ',mods{i}.mixer.desc]);
   end
end
disp(['npar = ',num2str(m.npar)])

m2 = Model2(frag,fnar,fdif);
pars = zeros(1,m2.npar);
m2.setPar(pars);
m2.sepKE = 1;
m2.sepSP = 0;
m2.rhodep = 0;

%%
clear classes;
load('temp.mat');
% KE diagonal    1: H  2: Cs  3: Cp
% Hen diagonal   4: H  5: Cs  6: Cp
% KE bonding     7: H-Cs  8: H-Cp  9: Cs-Cs 10: Cs-Cp 11: Cp-Cp  
% Hen bonding   12: H-Cs 13: H-Cp 14: Cs-Cs 15: Cs-Cp 16: Cp-Cp
m = Model3(frag,fnar,fdif);
m.addKEmodDiag(1,1);
m.addKEmodDiag(6,1);
m.addKEmodDiag(6,2);
m.addENmodDiag(1,1);
m.addENmodDiag(6,1);
m.addENmodDiag(6,2);

% m.addKEmodBonded(1,6,1,[1 2]);
% m.addKEmodBonded(6,6,[1 2],[1 2]);

 m.addKEmodBonded(1,6,1,1);
 m.addKEmodBonded(1,6,1,2);
 m.addKEmodBonded(6,6,1,1);
 m.addKEmodBonded(6,6,1,2);
 m.addKEmodBonded(6,6,2,2);

m.addENmodBonded(1,6,1,1);
m.addENmodBonded(1,6,1,2);
m.addENmodBonded(6,6,1,1);
m.addENmodBonded(6,6,1,2);
m.addENmodBonded(6,6,2,2);

% mods = m.ENmods;
% for i=1:size(mods,2)
%    disp(['KE ',num2str(i),' ilist [',num2str(mods{i}.ilist), ...
%       '] jlist [',num2str(mods{i}.jlist),'] ',mods{i}.mixer.desc]);
% end
for iatom = 1:m.natom
   disp(['H1 EN for atom ',num2str(iatom)]);
   mods = m.ENmods{1,iatom};
   for i=1:size(mods,2)
      disp(['EN ',num2str(i),' ilist [',num2str(mods{i}.ilist), ...
         '] jlist [',num2str(mods{i}.jlist),'] ',mods{i}.mixer.desc]);
   end
end

m2 = Model2(frag,fnar,fdif);
m2.sepKE = 1;
m2.sepSP = 1;
m2.rhodep = 0;

pars = rand(1,m2.npar+1);
m.setPar(pars);
m2.setPar(pars);
disp(['KE diff ',num2str( max(max(abs(m2.KE(1)-m.KE(1)))))]);

for iatom = 1:m.natom
   disp(['EN diff ',num2str(iatom),' ',num2str( max(max(abs(m2.H1en(iatom)-m.H1en(iatom)))))]);
end

%%
disp('m.solveHF');
m.solveHF;
%disp('m2.solveHF');
%m2.solveHF;
%disp(['orbital energy diffs ',num2str(max(max(abs(m.EorbEnv - m2.EorbEnv))))]);

%% Test of Model3.addENmodBonded1
clear classes;
load('temp.mat');
m1 = Model3(frag,fnar,fdif);
par1 = rand;
mix1 = Mixer(par1,1,'EN C H');
m1.addENmodBonded(1,6,1,[1 2],mix1);

m2 = Model3(frag,fnar,fdif);
mix2 = Mixer(par1,1,'EN1 C H');
mix3 = Mixer(par1,1,'EN1 H C');
m2.addENmodBonded1(1,6,1,[1 2],mix2);
m2.addENmodBonded1(6,1,[1 2],1,mix3);

enDiff = zeros(1,m1.natom);
for iatom = 1:m1.natom
   enDiff(iatom) = max(max(abs( m1.H1en(iatom,0) - m2.H1en(iatom,0))));
end

enDiff

%% Two electron modifications
clear classes;
% Need to have frag, fnar and fdif stored in temp.mat
load('temp.mat');
m = Model3(frag,fnar,fdif);
%m.addH2modDiag(1);
%m.addH2modDiag(6);
%m.addH2modOffDiag(1,1);
%m.addH2modOffDiag(1,6);
%m.addH2modOffDiag(6,6);
envs = 1:20;
E2f = frag.E2(envs);
m.solveHF(envs);
E2m = m.E2(envs);
E2diff = max(abs(E2f-E2m));

%% testing fitme
clear classes;
load('temp.mat');
m = Model3(frag,fnar,fdif);
mix1 = Mixer([0 0],4,'kec');
mix1.fixed = [0 1];
m.addH2modDiag(1);
m.addH2modDiag(6);
f1 = Fitme;
f1.addFrag(m,fhl);
f1.includeKE = 1;
f1.includeEN = zeros(1,6);
f1.includeE2 = 1;
f1.setEnvs(1:10);

m2 = m;
mfile = m;
save('scratch/junk1.mat','mfile');
load('scratch/junk1.mat');
m3 = mfile;


%%
start = f1.getPars;
limits = [];
options = optimset('DiffMinChange',1.0e-5,'TolFun',1.0e-4,'TolX',1.0e-3);
[pt,resnorm,residual,exitflag,output,lambda,jacobian] = ...
   lsqnonlin(@f1.err, start,-limits,limits,options);
