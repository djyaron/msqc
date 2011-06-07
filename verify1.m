%% Test of Fragment class
clear classes;
c1 = Fragment.defaultConfig();
c1.template = 'h2';
c1.basisSet = 'STO-3G';
i=0;
for bl=1.0:0.1:1.4
   i = i+1;
   disp(['bondlength = ', num2str(bl)]);
   c1.par = bl;
   frag(i,1) = Fragment('data',c1);
   H1t = sum(frag(i,1).H1en,3);
   r1 = H1t - (frag(i,1).H1 - frag(i,1).KE);
   disp(['sum of H1en agrees with H1 to ',num2str(max(max(abs(r1))))]);
end

%% Test of Environment class
clear classes;
size = [3,3,3];
e1 = Environment.newCube(size,1);
e1.plotFig(3);
e2 = Environment.newCube(size,1);

save('junk1.mat', 'e1','e2');

%% Do some calculations in environments
clear classes;
load('junk1.mat');
c1 = Fragment.defaultConfig();
c1.template = 'h2';
c1.basisSet = '6-31G**';
c1.par = 1.0;
frag = Fragment('data',c1);
%%
frag.setEnvSize(3);
frag.addEnv(e1);
frag.addEnv(e2);
frag.addEnv(e1);

%% test of load all environments
clear classes;
c1 = Fragment.defaultConfig();
c1.template = 'h2';
c1.basisSet = 'STO-3G';
c1.par = 1.1;
frag = Fragment('data',c1);
frag.loadAllEnv();

%% Test of Hartree fock routine
clear all;
c1 = Fragment.defaultConfig();
c1.template = 'h2';
c1.basisSet = '6-31G';
c1.par = 1.1;
frag = Fragment('data',c1);
frag.loadAllEnv();

for ienv = 0:frag.nenv
   [orb,Eorb,Ehf] = hartfock(frag,ienv);
   if (ienv == 0)
      Ecomp = frag.Eorb;
   else
      Ecomp = frag.EorbEnv(:,ienv);
   end
   Ediff = max(max(abs(Eorb - Ecomp)));
   disp(['env = ',num2str(ienv),' energy diff ',num2str(Ediff)]);
end

%% Test of partitioning
clear classes;
c1 = Fragment.defaultConfig();
c1.template = 'h2';
c1.basisSet = '6-31G';
c1.par = 1.1;
frag = Fragment('data',c1);
frag.loadAllEnv();

for ienv = 0:frag.nenv
   E2a = sum(sum(sum(sum(frag.density2p(ienv).*frag.H2))));
   if (ienv == 0)
      E1a = frag.Ehf - E2a - frag.Hnuc;
   else
      E1a = frag.EhfEnv(ienv) - E2a - frag.HnucEnv(ienv);
   end
   E1p = frag.partitionE1(ienv);
   E1t = sum(sum(E1p));
   disp(['env',num2str(ienv),':  E1 from gaussian ',num2str(E1a),' E1 from partition ', ...
      num2str(E1t),' diff ', num2str(E1a - E1t)]);
   E2p = frag.partitionE2(ienv);
   E2t = sum(sum(sum(sum(E2p))));
   disp(['env',num2str(ienv),':  E2 from gaussian ',num2str(E2a),' E2 from partition ', ...
   num2str(E2t),' diff ', num2str(E2a - E2t)]);
end

%% Test of Model1
clear classes;
c1 = Fragment.defaultConfig();
c1.template = 'h2';
c1.basisSet = '6-31G';
c1.par = 1.1;
frag = Fragment('data',c1);
frag.loadAllEnv();

m1 = Model1(frag);
m1.solveHF;

disp(['HF energies of isolated fragment agree to: ',...
   num2str(abs(frag.Ehf-m1.Ehf))]);
disp(['HF energies in environments agree to: ', ...
   num2str(max(abs(frag.EhfEnv - m1.EhfEnv)))]);
disp(['HF orbital energies agree to: ',...
   num2str(max(max(abs(frag.EorbEnv - m1.EorbEnv))))]);


