%% CH4 charge plots, comparison of high-level and low-level
clear all;
load('datasets/ethanerDat.mat');
qs = cell(6,2);  % {Z, LL or HL}
ch3sum = cell(2,1);
close all;
col = {'b','r'};
egs = cell(2,1);
for ilevel = 1:2
   for imod = 1:size(LL,1);
      switch ilevel
         case 1
            m1 = LL{imod,1};
         case 2
            m1 = HL{imod,1};
      end
      for ienv = 1:m1.nenv
         t1 = m1.mcharge(ienv)-m1.mcharge(0);
         egs{ilevel}(end+1) = m1.EhfEnv(ienv) - m1.Ehf;
         ch3sum{ilevel}(end+1) = sum(t1([1 3 4 5]));
         ch3sum{ilevel}(end+1) = sum(t1([2 6 7 8]));
         for iatom = 1:length(t1)
            qs{m1.Z(iatom),ilevel}(end+1) = t1(iatom);
         end
      end
   end
   figure(101)
   hold on;
   subplot(2,2,2*(ilevel-1) +1);
   hist(qs{1,ilevel},10);
   hold on;
   subplot(2,2,2*(ilevel-1) + 2);
   hist(qs{6,ilevel},10);
end
figure(300)
hist(ch3sum{1},20)
figure(301)
hist(ch3sum{2},20)
%%
disp(['mean of C charges from LL and HL ', num2str(mean(qs{6,1})), ' ', num2str(mean(qs{6,2}))]);
disp(['mean of H charges from LL and HL ', num2str(mean(qs{1,1})), ' ', num2str(mean(qs{1,2}))]);
disp(['std of C charges from LL and HL  ', num2str(std(qs{6,1})), ' ', num2str(std(qs{6,2}))]);
disp(['std of H charges from LL and HL  ', num2str(std(qs{1,1})), ' ', num2str(std(qs{1,2}))]);
disp(['std of CH3 charges from LL and HL  ', num2str(std(ch3sum{1})), ' ', num2str(std(ch3sum{2}))]);

%%  Effect of a point charge on methane
% make environments
envs = cell(1,1);
ic = 0;
% for r=2:0.25:5
%    q = 1;
%    ic = ic+1;
%    envs{ic} =Environment;
%    envs{ic}.addCharge([0 0 r], q);
% end
r = 5;
for q=1:0.5:4
   ic = ic+1;
   envs{ic} =Environment;
   envs{ic}.addCharge([0 0 r], q);
end
dataroot = 'c:\matdl\yaron\10-30-12\methane-q1\';
if (exist(dataroot,'dir') ~= 7)
   status = mkdir(dataroot);
end
copyfile('templates/ch4.tpl',[dataroot,'/ch4.tpl']);

config = Fragment.defaultConfig();
config.template = 'ch4';
config.par = [1.12 1.12 1.12 1.12 109.47 109.47 109.47 120.0 -120.0]


config.basisSet = 'STO-3G';
LL = Fragment([dataroot], config);
for ienv = 1:length(envs)
   display(['LL env ',num2str(ienv)]);
   LL.addEnv(envs{ienv});
end
config.basisSet = '6-31G';
HL = Fragment([dataroot], config);
for ienv = 1:length(envs)
   display(['HL env ',num2str(ienv)]);
   HL.addEnv(envs{ienv});
end
%%
close all
for i= 1:length(envs)
   %x = envs{i}.r(3,1);
   x = envs{i}.rho(1);
   figure(500);
   hold on;
   plot(x,(LL.EhfEnv(i) - LL.Ehf) * 627.509,'bo');
   plot(x,(HL.EhfEnv(i) - HL.Ehf) * 627.509,'ro');
   figure(600);
   hold on;
   qLL = LL.mcharge(i) - LL.mcharge(0);
   qHL = HL.mcharge(i) - HL.mcharge(0);
   plot(x,qLL(1),'bo');
   plot(x,qHL(1),'ro');
   plot(x,qLL(2),'bx');
   plot(x,qHL(2),'rx');
   plot(x,qLL(3),'b+');
   plot(x,qHL(3),'r+');
end

%% Effects of dipoles from cube with edges of 5
mags = [];
edge = 5;
ic = 0;
envs = cell(1,1);
for dip = 5:5:20
   for iavg = 1:10
      ic = ic + 1;
      envs{ic} = Environment.dipCube(edge,0.1,dip);
      mags(ic) = dip;
   end
   dataroot = 'c:\matdl\yaron\10-30-12\methane-dip';
   if (exist(dataroot,'dir') ~= 7)
      status = mkdir(dataroot);
   end
end
save([dataroot,'\envs.mat'], 'envs','mags')
copyfile('templates/ch4.tpl',[dataroot,'\ch4.tpl']);
%%
clear all
dataroot = 'c:\matdl\yaron\10-30-12\methane-dip';
load([dataroot,'\envs.mat'], 'envs','mags');
config = Fragment.defaultConfig();
config.template = 'ch4';
config.par = [1.12 1.12 1.12 1.12 109.47 109.47 109.47 120.0 -120.0];

config.basisSet = 'STO-3G';
LL = Fragment([dataroot], config);
for ienv = 1:length(envs)
   display(['LL env ',num2str(ienv)]);
   LL.addEnv(envs{ienv});
end
config.basisSet = '6-31G';
HL = Fragment([dataroot], config);
for ienv = 1:length(envs)
   display(['HL env ',num2str(ienv)]);
   HL.addEnv(envs{ienv});
end
%%
close all
xold = mags(1);
yLL = [];
yHL = [];
for i = 1:length(envs);
   x = mags(i);
   figure(500);
   hold on;
   plot(x,(LL.EhfEnv(i) - LL.Ehf) * 627.509,'bo');
   plot(x,(HL.EhfEnv(i) - HL.Ehf) * 627.509,'ro');
   figure(600);
   hold on;
   qLL = LL.mcharge(i) - LL.mcharge(0);
   qHL = HL.mcharge(i) - HL.mcharge(0);
   if ((x == xold) && i<length(envs))
      yLL(end+1) = qLL(1);
      yHL(end+1) = qHL(1);
   else
      plot(xold,std(yLL),'bo');
      plot(xold,std(yHL),'ro');
      xold = x; yLL = []; yHL = [];
   end
%    plot(x,qLL(2),'bx');
%    plot(x,qHL(2),'rx');
%    plot(x,qLL(3),'b+');
%    plot(x,qHL(3),'r+');
end
%% Generate a new set of environments
clear all;
nenv = 25;
edge = 5;
env = cell(25,1);
for ienv = 1:25
   t1 = Environment.dipCube(edge,0.1,25);
   r = -1 + (2 *rand(3,1));  % each element will go from -1 to +1 
   orien = r/norm(r);
   orien = 5 * orien;
   rho = 5 * (-1 + 2 * rand(1,1));
   t1.addCharge(orien,rho);
   env{ienv} = t1;
end
save('datasets/env3.mat','env');
%% Generate a new set of environments
clear all;
nenv = 25;
edge = 5;
rhos = -5:0.5:5;
env = cell(length(rhos),1);
ienv = 0;
for rho=rhos
   ienv = ienv+1;
   t1 = Environment.dipCube(edge,0.1,25);
   r = -1 + (2 *rand(3,1));  % each element will go from -1 to +1 
   orien = [0 0 5];
   t1.addCharge(orien,rho);
   env{ienv} = t1;
end
save('datasets/env4.mat','env');
%% Generate a new set of environments
clear all;
edge = 5;
env = cell(25,1);
for ienv = 1:25
   env{ienv} = Environment.dipCube(edge,0.1,25);
end
save('datasets/env6.mat','env');
%% Generate a new set of environments for ethane
clear all;
nenv = 25;
ccBond = 1.54; % bond length
edge = [5, 5, 5+(ccBond/2)];
env = cell(nenv,1);
for ienv = 1:nenv
   t1 = Environment.dipCube(edge,0.1,25);
   r = -1 + (2 *rand(3,1));  % each element will go from -1 to +1 
   orien = [0 0 (5+(ccBond/2))];
   rho = 5 * (-1 + 2 * rand(1,1));
   t1.addCharge(orien,rho);
   env{ienv} = t1;
end
save('datasets/ethane-env.mat','env');


%% Look at this data
clear;
load('datasets/env3.mat','env');
load('datasets/ch4rDat.mat');
%%
rho = [];
for ienv = 1:length(env)
   rho(ienv) = env{ienv}.rho(end);
end
[rs,is] = sort(rho);
%%
close all;
for ienv = 1:length(env)
   hold on;
   env{ienv}.plotFig(1);
   junk = input('hit next');
end
%%
qs = cell(6,2);  % {Z, LL or HL}
ch3sum = cell(2,1);
close all;
col = {'b','r'};
lev = {'LL','HL'};
egs = cell(2,1);
qc = [];
for ilevel = 1:2
   for imod = 1:size(LL,1);
      switch ilevel
         case 1
            m1 = LL{imod,1};
         case 2
            m1 = HL{imod,1};
      end
      for ienv = 1:m1.nenv
         t1 = m1.mcharge(ienv)-m1.mcharge(0);
         egs{ilevel}(end+1) = m1.EhfEnv(ienv) - m1.Ehf;
         %ch3sum{ilevel}(end+1) = sum(t1([1 3 4 5]));
         %ch3sum{ilevel}(end+1) = sum(t1([2 6 7 8]));
         qc(imod,ienv,ilevel) = t1(1);
         for iatom = 1:length(t1)
            qs{m1.Z(iatom),ilevel}(end+1) = t1(iatom);
         end
      end
   end
   figure(101)
   hold on;
   subplot(2,2,2*(ilevel-1) +1);
   hist(qs{1,ilevel},10);
   title(['H  ',lev{ilevel}]);
   hold on;
   subplot(2,2,2*(ilevel-1) + 2);
   hist(qs{6,ilevel},10);
   title(['C  ',lev{ilevel}]);
end
% figure(300)
% hist(ch3sum{1},20)
% figure(301)
% hist(ch3sum{2},20)
%%
close all;
for imod = 1:size(LL,1)
   t1 = qc(imod,:,1);
   t2 = qc(imod,:,2);
   figure(600)
   hold on;
   plot(t1,t2,'r.');
   figure(601)
   hold on;
   plot(rho,abs(t2-t1),'b.');
   figure(602)
   hold on;
   figure(603)
   hold on;
   plot(rho,t1,'b.');
   plot(rho,t2,'r.');
end
%%
disp(['mean of C charges from LL and HL ', num2str(mean(qs{6,1})), ' ', num2str(mean(qs{6,2}))]);
disp(['mean of H charges from LL and HL ', num2str(mean(qs{1,1})), ' ', num2str(mean(qs{1,2}))]);
disp(['std of C charges from LL and HL  ', num2str(std(qs{6,1})), ' ', num2str(std(qs{6,2}))]);
disp(['std of H charges from LL and HL  ', num2str(std(qs{1,1})), ' ', num2str(std(qs{1,2}))]);
%disp(['std of CH3 charges from LL and HL  ', num2str(std(ch3sum{1})), ' ', num2str(std(ch3sum{2}))]);
 

%% What is our current spread in energies
load('datasets\ch4rDat-orig.mat');
%load('datasets\ethylenerDat.mat');
%%
nmol = size(LL,1);
nenv = LL{1,1}.nenv; 
% +1 for ground state
Emol = zeros(nmol,nenv+1);
Eenv = Emol;
Egaussian = Emol;
for i=1:nmol
   f1 = LL{i,1};
   Eke1 = f1.EKE;
   Ehf = Eke1;
   for iatom = 1:f1.natom
      Een1{iatom} = f1.Een(iatom);
      Ehf = Ehf + Een1{iatom};
   end
   for ienv=0:f1.nenv;
      E21(ienv+1) = f1.partitionE2(ienv,f1.H2,{1:f1.nbasis});
   end
   Ehf = Ehf + E21;
   Ehf = Ehf + f1.Hnuc;
   Emol(i,:) = Ehf;
   Eenv(i,:) = f1.Eenv;
   % The [] makes a vector with first element no env, and remainder in env
   Egaussian(i,:) = [f1.Ehf, f1.EhfEnv];
end
Emol = Emol * 627.509;
Eenv = Eenv * 627.509;
Egaussian = Egaussian * 627.509;

% Emol is energy of molecule (without interaction with env) (geom,env)

disp(['max disagreement is ', ...
   num2str(max(max(abs(Egaussian - Emol - Eenv))))]);
figure(1);
hist(Emol(:,1)-mean(Emol(:,1)),20);
title('molecule in zero environment');
figure(2);
hist(Emol(:)-mean(Emol(:)),20);
title('molecule in environments');
Epol = zeros(size(Emol));
for i=1:size(Emol,1)
   Epol(i,:) = Emol(i,:) - Emol(i,1);
end
figure(3);
hist(Epol(:), 50);
title('polarization effects');
figure(8);
boxplot(Epol);
%%
figure(9);
boxplot(Epol');
